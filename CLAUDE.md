# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with the EntrapmentAnalyses.jl module.

## Project Overview

EntrapmentAnalyses.jl is a Julia package for entrapment-based false discovery rate (FDR) analysis in proteomics data. It implements both combined and paired empirical FDR methods, supporting both precursor-level and protein-level analyses.

**Key Features:**
- Parquet and TSV data loading
- Combined and paired EFDR calculations
- Per-file q-value calculation
- Protein-level rollup and analysis
- Notebook-style visualization with specific color schemes
- Progress monitoring for O(n²) operations

## Important Implementation Details

### 1. Critical Bug Fix in Paired EFDR
The paired EFDR calculation was fixed to match the notebook's exact logic using mutually exclusive conditions:
```julia
# Correct implementation with elseif:
if e_score >= s && s > o_score
    Nεsτ += 1
elseif e_score > o_score && o_score >= s
    Nετs += 1
end
```

### 2. Per-File Q-value Calculation
**IMPORTANT**: Q-values must be calculated separately for each file, not on concatenated data:
- Group by `file_name` first
- Calculate q-values within each group
- Use these per-file q-values for EFDR calculations

### 3. Plex-Specific Pairing Logic
**IMPORTANT**: The pairing system now matches the notebook's behavior:
- Complement scores are calculated per file AND per plex
- The same peptide in different plexes can have different complement scores
- Dictionary key is `(plex, pair_id)` not just `pair_id`
- Process: Build global dictionaries from library → Process each file separately → Build plex score dictionary per file

### 4. Data Format Expectations
**PSM Results (Parquet):**
- `stripped_seq`: Peptide sequence (String)
- `z`: Charge state (numeric)
- `PredVal`: Prediction score (Float)
- `decoy`: Decoy status (Bool)
- `file_name`: Source file identifier (String)
- `protein`: Protein identifiers (String)

**Spectral Library (TSV):**
- `PeptideSequence`: Modified peptide sequence (String)
- `PrecursorCharge`: Charge state (Int)
- `EntrapmentGroupId`: 0 for originals, >0 for entrapments (Int)
- `PrecursorIdx`: Unique pair identifier (Int)

## Common Development Commands

### Julia Package Management
```bash
# Enter Julia REPL
julia

# Activate the project environment
] activate .

# Install dependencies
] instantiate

# Add a new dependency
] add PackageName
```

### Running Tests
```bash
# From Julia REPL with activated environment
] test

# Or run specific test files
julia --project=. test/unit/test_efdr_methods.jl
```

### Development Workflow
```bash
# Start Julia with Revise for auto-reloading
julia --project=.

# In the REPL
using Revise, EntrapmentAnalyses

# Load example data
psm_df = load_parquet("path/to/psm_results.parquet")
spec_lib_df = load_spectral_library("path/to/spectral_library.tsv")

# Run analysis
results = run_efdr_analysis([psm_file], library_file; output_dir="output")
```

## Main API Functions

### 1. `run_efdr_analysis`
```julia
run_efdr_analysis(parquet_files::Vector{String}, library_path::String; kwargs...)
```
- Runs precursor-level EFDR analysis
- Calculates both combined and paired EFDR by default
- Outputs TSV results and plots

### 2. `run_protein_efdr_analysis`
```julia
run_protein_efdr_analysis(parquet_files::Vector{String}, library_path::String; kwargs...)
```
- Runs protein-level EFDR analysis
- Performs protein rollup (best PSM per protein)
- Calculates per-run protein q-values

## Code Architecture

```
src/
├── EntrapmentAnalyses.jl      # Main module, exports
├── io/
│   └── data_loaders.jl        # Parquet/TSV loading with missing value handling
├── core/
│   ├── pairing.jl             # Plex-aware pairing system with type definitions
│   ├── efdr_methods.jl        # EFDR calculations (unchanged)
│   └── scoring.jl             # Monotonization
├── analysis/
│   ├── qvalue_calculation.jl  # Q-value calculations
│   ├── protein_analysis.jl    # Protein rollup
│   └── efdr_analysis.jl       # Analysis coordination
├── plotting/
│   └── visualization.jl        # Notebook-style plots
└── api.jl                     # Public API functions (updated for plex pairing)
```

## Testing Guidelines

### Key Test Areas
1. **Paired EFDR Logic**: Test the fixed mutually exclusive conditions
2. **Per-File Q-values**: Ensure q-values are calculated per file
3. **Protein Rollup**: Verify best PSM selection per protein
4. **Plot Styling**: Check exact RGB values and sizes

### Running Specific Tests
```julia
# Test the critical paired EFDR fix
include("test/unit/test_paired_efdr_validation.jl")

# Test protein analysis
include("test/unit/test_protein_analysis.jl")
```

## Key Implementation Functions

### Pairing System
1. **`init_entrapment_pairs_dict`**: Builds global dictionaries from library
   - Returns `pair_dict` and `is_original_dict`
   - Called once per analysis

2. **`add_plex_complement_scores!`**: Adds plex-specific complement scores
   - Processes each file separately
   - Builds `plex_prec_to_scores` dictionary per file
   - Adds `complement_score`, `is_original`, and `pair_id` columns

3. **`compute_pairing_vectors!`**: Main pairing function (mutating)
   - Calls the above functions
   - Modifies DataFrame directly by adding columns
   - No longer returns pairing_info tuple

### Data Loading
1. **`handle_missing_values!`**: Centralized missing value handling
   - Shows accurate warnings with column names and counts
   - Validates types after replacement

## Common Issues and Solutions

### Issue: "UndefVarError: load_parquet not defined"
**Solution**: The function is exported, but Julia needs to be restarted after module changes.

### Issue: Q-values not monotonic
**Solution**: Ensure `monotonize!` is called after q-value calculation and operates on sorted data.

### Issue: Missing complement pairs
**Solution**: The pairing system handles unpaired sequences by setting complement_indices to -1.

### Issue: Type instability warnings
**Solution**: Use type assertions when accessing DataFrame columns:
```julia
scores = df[!, :PredVal]::AbstractVector{Float32}
```

### Issue: "TypeError: expected Vector{Float32}, got SubArray"
**Solution**: When working with grouped DataFrames, use `AbstractVector` instead of `Vector` for type assertions:
```julia
# Wrong - causes TypeError with SubDataFrames
score_vec = group[!, score_col]::Vector{Float32}

# Correct - handles both Vector and SubArray types
score_vec = group[!, score_col]::AbstractVector{Float32}
```

## Visualization Specifications

**Exact Plot Styling (from notebook):**
- Size: 600×450 pixels (`size = (400*1.5, 300*1.5)`)
- Color: RGB(0.39215686, 0.58431373, 0.92941176)
- Alpha: 0.75
- Line width: 3
- Font sizes: title=16, axis labels=16, ticks=12
- Axis limits: (0, 0.05) for both x and y

## Performance Considerations

1. **Paired EFDR is O(n²)**: Progress bars show actual operation count
2. **Pre-computed pairing**: Avoids dictionary lookups in main loops
3. **Vector operations**: Use typed vectors instead of DataFrame columns in loops
4. **Type stability**: Ensure all function arguments have known types

## Recent Changes

### Q-value Calculation Fix
- Changed from calculating q-values on concatenated data to per-file calculation
- Added `calculate_qvalues_per_file!` function in `qvalue_calculation.jl`
- Fixed MethodError by changing function signatures to accept `AbstractDataFrame`

### Performance Enhancement
- Added local q-value filtering (default 0.05) to reduce computation time
- Added `local_qval_threshold` parameter to `run_efdr_analysis`

### Visualization Updates
- X-axis limited to `min(local_qval_threshold, global_qval_threshold)`
- Y-axis dynamically scaled based on data
- Updated all plotting functions to accept `xlim` parameter

### EFDR Calculation Fixes
- Fixed paired EFDR assignment issue - now assigns directly to parent DataFrame indices
- Fixed combined EFDR to be calculated per-run instead of globally
- Both combined and paired EFDR are now calculated inside the same per-file loop
- Added debug output for both EFDR methods showing max values per run

### Report Generation
- Added comprehensive markdown report generation with `generate_analysis_report`
- Generates separate plots for combined EFDR, paired EFDR, and comparison
- Saves plots as both PDF and PNG for markdown embedding

### Type Assertion Fix
- Changed all `Vector{T}` type assertions to `AbstractVector{T}`
- Fixes TypeError when working with SubDataFrames from grouped operations
- Allows handling both Vector and SubArray types from DataFrame column access

### Plex-Specific Pairing Implementation
- Implemented plex-aware pairing system that matches notebook behavior
- Added new type definitions: `PeptideKey`, `PlexPairKey`, `ScorePair`
- Created `compute_pairing_vectors` that returns plex-specific complement scores
- Complement scores now respect both file and plex boundaries
- Deprecated `get_complement_scores` in favor of pre-computed plex-specific values

### Data Loader Improvements
- Fixed incorrect warning messages in `load_parquet_results` and `load_spectral_library`
- Added `handle_missing_values!` helper function to eliminate code duplication
- Defined column specifications with proper types and descriptions
- All columns now properly checked for missing values with accurate warnings
- Shows count of missing values for better debugging

### Pairing API Refactor
- Changed `compute_pairing_vectors` to `compute_pairing_vectors!` (mutating function)
- Function now modifies DataFrame directly by adding columns:
  - `is_original`: Boolean flag (true=original, false=entrapment)
  - `pair_id`: Links original/entrapment pairs from the library
  - `entrap_label`: Integer label for EFDR calculation (0=original, 1=entrapment)
  - `complement_score`: Plex-specific score of the paired peptide
  - `complement_indices`: Deprecated, kept for compatibility
- Removed `pairing_info` named tuple return value
- Updated all API functions to use DataFrame columns directly
- This simplifies the API and makes data flow more transparent

## Development Guidelines for Claude

### Key Instructions
1. **Make commits as you go** - Commit changes incrementally as you work on tasks
2. **Don't run tests yourself** - Ask the user to run tests and provide instructions on how to run them
3. **Ask for help when confused** - Request clarification whenever there's uncertainty about requirements or implementation

## Future Improvements to Consider

1. Parallelize per-file EFDR calculations
2. Add streaming support for very large datasets
3. Implement combined+paired plot overlays
4. Add statistical summaries to output reports