# EntrapmentAnalyses.jl

A Julia module for entrapment-based false discovery rate (FDR) analysis in proteomics data.

## Installation

```julia
using Pkg
Pkg.add(path="/path/to/EntrapmentAnalyses")
```

## Quick Start

```julia
using EntrapmentAnalyses

# Example file paths
parquet_file = "/Users/nathanwamsley/Data/May2025/kmd_jmod_search/9plex/05202025/all_IDs_filtered_01.parquet"
library_file = "/Users/nathanwamsley/Data/May2025/spec_libs/parsed_libs/hs_tag6_predlib_JDRT_480_1000_2ng_shufentrap_noloss_051925_jmod.tsv"

# Run precursor-level EFDR analysis (calculates both combined and paired by default)
results = run_efdr_analysis([parquet_file], library_file; output_dir="efdr_output")

# Run protein-level EFDR analysis
protein_results = run_protein_efdr_analysis([parquet_file], library_file; output_dir="protein_efdr_output")

# Generate comprehensive report with all visualizations
generate_analysis_report(results, "analysis_output")
```

## Key Functions

### Data Loading
- `load_parquet(filepath)` - Load PSM results from Parquet file
- `load_spectral_library(filepath)` - Load spectral library from TSV file

### Analysis Functions
- `run_efdr_analysis(parquet_files, library_path)` - Precursor-level EFDR analysis
- `run_protein_efdr_analysis(parquet_files, library_path)` - Protein-level rollup and analysis

### Visualization
- `plot_combined_efdr(df)` - Plot combined EFDR results
- `plot_paired_efdr(df)` - Plot paired EFDR results
- `plot_efdr_comparison_both_methods(df)` - Compare both EFDR methods
- `generate_analysis_report(df, output_dir)` - Generate comprehensive markdown report

## Data Format

### PSM Results (Parquet)
Required columns:
- `SpecId`: Spectrum identifier
- `Label`: Multiplex channel label
- `Peptide`: Peptide sequence
- `Proteins`: Protein identifiers
- `PredVal`: Prediction score
- `ExpMass`: Experimental mass
- `CalcMass`: Calculated mass

### Spectral Library (TSV)
Required columns:
- `SpecId`: Spectrum identifier
- `ModifiedPeptide`: Modified peptide sequence
- `PrecursorCharge`: Charge state
- `iRT`: Retention time
- `PrecursorMz`: Precursor m/z

## Advanced Usage

### Working with Multiple Files

```julia
# Process multiple Parquet files
parquet_files = [
    "file1.parquet",
    "file2.parquet",
    "file3.parquet"
]
results = run_efdr_analysis(parquet_files, library_file)
```

### Custom Parameters

```julia
# Run with custom parameters
results = run_efdr_analysis(
    parquet_files, 
    library_file;
    output_dir = "custom_output",
    global_qval_threshold = 0.05,  # Filter at 5% global FDR
    r_lib = 2.0,                   # Custom ratio parameter
    show_progress = false          # Disable progress bars
)

# Generate report with custom column names
generate_analysis_report(
    results,
    "report_output";
    combined_efdr_col = :combined_entrapment_fdr,
    paired_efdr_col = :precursor_entrapment_fdr
)
```

## Algorithm Details

The module implements entrapment-based FDR estimation following the Noble lab methodology:

- **Combined EFDR**: `(Nε × (1 + 1/r)) / (Nτ + Nε)`
- **Paired EFDR**: `(Nε + Nεsτ + 2×Nετs) / (Nτ + Nε)`

Where:
- `Nτ`: Number of targets above threshold
- `Nε`: Number of entrapments above threshold  
- `Nεsτ`: Entrapments winning with score ≥ threshold > paired target
- `Nετs`: Entrapments winning with both ≥ threshold