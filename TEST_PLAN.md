# EntrapmentAnalyses.jl Testing Plan

## Overview
This document outlines the comprehensive testing plan for the EntrapmentAnalyses.jl package. Tests will be organized to mirror the src/ directory structure and will be executed using Julia's built-in package testing system (`]test` or `Pkg.test()`).

## Progress Tracking
- [ ] Step 1: Create dummy test and verify `]test` works
- [ ] Step 2: Set up test directory structure
- [ ] Step 3: Test io/data_loaders.jl
- [ ] Step 4: Test core/pairing.jl
- [ ] Step 5: Test core/efdr_methods.jl
- [ ] Step 6: Test core/scoring.jl
- [ ] Step 7: Test analysis/qvalue_calculation.jl
- [ ] Step 8: Test analysis/protein_analysis.jl
- [ ] Step 9: Test analysis/efdr_analysis.jl
- [ ] Step 10: Test plotting/visualization.jl
- [ ] Step 11: Integration tests for api.jl
- [ ] Step 12: Final test suite validation

## Test Structure
```
test/
├── runtests.jl              # Main test runner
├── unit/                    # Unit tests
│   ├── io/
│   │   └── test_data_loaders.jl
│   ├── core/
│   │   ├── test_pairing.jl
│   │   ├── test_efdr_methods.jl
│   │   └── test_scoring.jl
│   ├── analysis/
│   │   ├── test_qvalue_calculation.jl
│   │   ├── test_protein_analysis.jl
│   │   └── test_efdr_analysis.jl
│   └── plotting/
│       └── test_visualization.jl
├── integration/             # Integration tests
│   └── test_api.jl
└── test_data/              # Test fixtures
    ├── sample_psm.parquet
    ├── sample_library.tsv
    └── expected_outputs/
```

## Step 1: Dummy Test Setup

### Objective
Create a minimal test that verifies the Julia testing infrastructure works correctly.

### Implementation
1. Create `test/runtests.jl` with a simple test
2. Verify package can be loaded
3. Run a basic assertion

### Files to Create
- `test/runtests.jl`

## Step 2: Test Directory Structure

### Objective
Set up the complete test directory structure mirroring src/.

### Implementation
1. Create all necessary subdirectories
2. Create placeholder test files
3. Update runtests.jl to include all test files

## Step 3: Test io/data_loaders.jl

### Functions to Test
- `load_parquet_results()`
- `load_spectral_library()`
- `handle_missing_values!()`

### Test Cases
1. **load_parquet_results**
   - Loading valid parquet file
   - Handling missing values in each column
   - Type conversions (Int64 → Float32)
   - Empty file handling
   - Non-existent file error

2. **load_spectral_library**
   - Loading valid TSV file
   - Column renaming
   - Missing value handling
   - Invalid file format error

3. **handle_missing_values!**
   - Replace missing values correctly
   - Type validation after replacement
   - Warning message accuracy

## Step 4: Test core/pairing.jl

### Functions to Test
- `init_entrapment_pairs_dict()`
- `add_plex_complement_scores!()`
- `compute_pairing_vectors!()`

### Test Cases
1. **init_entrapment_pairs_dict**
   - Correct dictionary construction
   - Original vs entrapment identification
   - Handling unpaired peptides

2. **add_plex_complement_scores!**
   - Plex-specific score assignment
   - Handling missing pairs
   - Column addition to DataFrame

3. **compute_pairing_vectors!**
   - Full pipeline test
   - Column presence validation
   - Data integrity checks

## Step 5: Test core/efdr_methods.jl

### Functions to Test
- `calculate_combined_efdr!()`
- `calculate_paired_efdr!()`

### Test Cases
1. **calculate_combined_efdr!**
   - Correct EFDR calculation
   - Handling edge cases (no entrapments, all entrapments)
   - Per-file calculation

2. **calculate_paired_efdr!**
   - Mutually exclusive condition logic
   - O(n²) performance with small datasets
   - Correct EFDR assignment

## Step 6: Test core/scoring.jl

### Functions to Test
- `monotonize!()`

### Test Cases
1. **monotonize!**
   - Non-increasing scores maintained
   - Already monotonic data unchanged
   - Empty vector handling
   - Single element vector

## Step 7: Test analysis/qvalue_calculation.jl

### Functions to Test
- `adjust_pvalues_bh()`
- `calculate_qvalues!()`
- `calculate_qvalues_per_file!()`

### Test Cases
1. **adjust_pvalues_bh**
   - Benjamini-Hochberg correctness
   - Monotonicity of output
   - Edge cases (all same p-value)

2. **calculate_qvalues!**
   - Q-value calculation accuracy
   - Monotonization applied
   - Column addition

3. **calculate_qvalues_per_file!**
   - Per-file grouping
   - Consistent with single-file calculation

## Step 8: Test analysis/protein_analysis.jl

### Functions to Test
- `get_best_psm_per_protein()`
- `calculate_protein_qvalues!()`

### Test Cases
1. **get_best_psm_per_protein**
   - Best PSM selection
   - Multiple proteins per PSM
   - Score handling

2. **calculate_protein_qvalues!**
   - Per-run q-value calculation
   - Protein-level FDR accuracy

## Step 9: Test analysis/efdr_analysis.jl

### Functions to Test
- `prepare_efdr_data()`
- `run_precursor_efdr_analysis()`
- `run_protein_efdr_analysis()`

### Test Cases
1. Integration of all components
2. Output file generation
3. Report generation

## Step 10: Test plotting/visualization.jl

### Functions to Test
- `plot_efdr_curves()`
- `save_plots()`

### Test Cases
1. Plot generation without errors
2. Correct styling (colors, sizes)
3. File output verification

## Step 11: Integration Tests for api.jl

### Functions to Test
- `run_efdr_analysis()`
- `run_protein_efdr_analysis()`

### Test Cases
1. Full pipeline execution
2. Single vs multi-file handling
3. All output files generated
4. Parameter validation

## Step 12: Final Validation

### Objectives
1. Run full test suite
2. Check code coverage
3. Performance benchmarks
4. Documentation completeness

## Running Tests

### Instructions for Running Tests

1. **Run all tests:**
   ```julia
   # From Julia REPL
   ] test
   
   # Or from Julia code
   using Pkg
   Pkg.test()
   ```

2. **Run specific test file:**
   ```julia
   # From project root
   julia --project=. test/unit/io/test_data_loaders.jl
   ```

3. **Run with coverage:**
   ```julia
   using Pkg
   Pkg.test(coverage=true)
   ```

## Test Data Requirements

### Sample Files Needed
1. `sample_psm.parquet` - Small PSM results file with known values
2. `sample_library.tsv` - Corresponding spectral library
3. Expected output files for validation

### Data Characteristics
- Include both original and entrapment peptides
- Multiple files for multi-file testing
- Edge cases (empty, single peptide, etc.)

## Notes
- After each step implementation, ask user to run `]test` and provide feedback
- Make commits after each successful test implementation
- Update this checklist as progress is made