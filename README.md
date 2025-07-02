# EntrapmentAnalyses.jl

[![Build Status](https://github.com/nathanwamsley/EntrapmentAnalyses.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/nathanwamsley/EntrapmentAnalyses.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![codecov](https://codecov.io/gh/nathanwamsley/EntrapmentAnalyses.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/nathanwamsley/EntrapmentAnalyses.jl)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://nathanwamsley.github.io/EntrapmentAnalyses.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://nathanwamsley.github.io/EntrapmentAnalyses.jl/dev)

A Julia package for calculating entrapment-based false discovery rate (FDR) in proteomics data. Implements both combined and paired empirical FDR methods for precursor-level and protein-level analyses. 

## Features

- **Dual EFDR Methods**: Implements both combined and paired empirical FDR calculations
- **Multi-level Analysis**: Supports both precursor-level and protein-level analyses
- **Plex-Aware Pairing**: Handles multiplexed data with plex-specific complement scoring
- **Per-File Processing**: Calculates q-values separately for each input file
- **Comprehensive Visualization**: Generates publication-ready plots with exact notebook styling
- **Detailed Reporting**: Creates markdown reports with analysis summaries and embedded plots

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

# Single file analysis
parquet_file = "path/to/psm_results.parquet"
library_file = "path/to/spectral_library.tsv"

# Run precursor-level EFDR analysis (calculates both combined and paired EFDR)
results = run_efdr_analysis(parquet_file, library_file; output_dir="efdr_output")

# Run protein-level EFDR analysis
protein_results = run_protein_efdr_analysis(parquet_file, library_file; output_dir="protein_output")

# Multiple file analysis
parquet_files = ["file1.parquet", "file2.parquet", "file3.parquet"]
results = run_efdr_analysis(parquet_files, library_file; output_dir="multi_file_output")
```

## Key Functions

### Main Analysis Functions
- `run_efdr_analysis(parquet_files, library_path; kwargs...)` - Precursor-level EFDR analysis
- `run_protein_efdr_analysis(parquet_files, library_path; kwargs...)` - Protein-level rollup and analysis

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
- `protein`: Protein identifiers
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
- `precursor_entrapment_fdr` or `protein_group_entrapment_fdr`: Paired EFDR values

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
# - Protein rollup and analysis
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