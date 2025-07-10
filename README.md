# EntrapmentAnalyses.jl

[![Build Status](https://github.com/nwamsley1/EntrapmentAnalyses.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/nwamsley1/EntrapmentAnalyses.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![codecov](https://codecov.io/gh/nwamsley1/EntrapmentAnalyses.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/nwamsley1/EntrapmentAnalyses.jl)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://nwamsley1.github.io/EntrapmentAnalyses.jl/dev)

EntrapmentAnalyses.jl is a Julia package that implementsw that combined and paired empirical FDR (EFDR) methods from [Wen et al. 2025](https://pubmed.ncbi.nlm.nih.gov/40524023/)


## Installation
Navigate to the EntrapmentAnalyses directory. Use `]` to enter Pkg mode. 
```julia

(@v1.11) pkg> activate .
  Activating project at `~/Projects/EntrapmentAnalysesJmod/EntrapmentAnalyses`

julia> using Revise, EntrapmentAnalyses
Precompiling EntrapmentAnalyses...
  1 dependency successfully precompiled in 6 seconds. 273 already precompiled.
```

## Quick Start

```julia
using EntrapmentAnalyses

# Run precursor-level EFDR analysis
results = run_efdr_analysis(
    ["data/psm_results.parquet"],
    "data/spectral_library.tsv";
    output_dir="output"
)
```

## Key Functions

### Main Analysis Functions
- `run_efdr_analysis(parquet_files, library_path; kwargs...)` - Precursor-level EFDR analysis

### Data Loading
- `load_parquet_results(filepaths)` - Load and combine PSM results from Parquet files
- `load_spectral_library(filepath)` - Load spectral library from TSV file

### Visualization
- `plot_combined_efdr(df)` - Plot combined EFDR results
- `plot_paired_efdr(df)` - Plot paired EFDR results  
- `plot_efdr_comparison_both_methods(df)` - Compare both EFDR methods
- `generate_analysis_report(df, output_dir)` - Generate comprehensive markdown report with all plots

## Data Format

### PSM Results (Parquet)
Required columns:
- `stripped_seq`: Peptide sequence without modifications
- `z`: Charge state
- `PredVal`: Prediction score (higher is better)
- `decoy`: Boolean indicating decoy status
- `file_name`: Source file identifier
- `channel`: Multiplex channel (optional, will add dummy if missing)

### Spectral Library (TSV)
Required columns:
- `PeptideSequence`: Modified peptide sequence
- `PrecursorCharge`: Charge state
- `EntrapmentGroupId`: 0 for originals, >0 for entrapments
- `PrecursorIdx`: Unique pair identifier linking originals to entrapments

## Output Files

Each analysis run generates:
- **TSV Results**: Complete results with EFDR columns added
- **PDF/PNG Plots**: Visualization of FDR vs EFDR curves
- **Markdown Report**: Comprehensive analysis summary with embedded plots

### Output Columns Added
- `local_qvalue`: Per-file q-values
- `global_qvalue`: Global q-values (best per precursor)
- `is_original`: Boolean indicating if peptide is original
- `pair_id`: Links original/entrapment pairs
- `complement_score`: Plex-specific score of paired peptide
- `combined_entrapment_fdr`: Combined EFDR values
- `precursor_entrapment_fdr`: Paired EFDR values

## Advanced Usage

### Custom Parameters

```julia
# Run with filtering and custom parameters
results = run_efdr_analysis(
    parquet_files, 
    library_file;
    output_dir = "custom_output",
    global_qval_threshold = 0.05,     # Filter at 5% global FDR
    local_qval_threshold = 0.01,      # Filter at 1% local FDR
    r_lib = 2.0,                      # Library to real entrapment ratio
    show_progress = false             # Disable progress bars
)

# Generate report with custom x-axis limit
generate_analysis_report(
    results,
    "report_output";
    xlim = (0, 0.1)  # Show FDR up to 10%
)
```

## Testing

The package includes a comprehensive test suite covering all major functionality:

```julia
# Run all tests
] test

# Tests cover:
# - Data loading with missing value handling
# - Plex-aware pairing system
# - Combined and paired EFDR calculations
# - Q-value calculations (per-file and global)
# - Monotonization functions
# - Visualization functions
```

## Algorithm Details

The module implements entrapment-based FDR estimation:

### Combined EFDR
```
EFDR = (Nε × (1 + 1/r)) / (Nτ + Nε)
```

### Paired EFDR
```
EFDR = (Nε + Nεsτ + 2×Nετs) / (Nτ + Nε)
```

Where:
- `Nτ`: Number of targets above threshold
- `Nε`: Number of entrapments above threshold  
- `Nεsτ`: Entrapments winning with score ≥ threshold > paired target
- `Nετs`: Entrapments winning with both ≥ threshold

## Requirements

- Julia 1.10+
- Dependencies: DataFrames, CSV, Parquet, ProgressBars, Plots, Dates

## License

MIT License