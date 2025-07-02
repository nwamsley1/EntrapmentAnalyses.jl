# User Guide

This guide provides comprehensive instructions for using EntrapmentAnalyses.jl in your proteomics workflows.

## Installation

EntrapmentAnalyses.jl requires Julia 1.0 or later. Navigate to the EntrapmentAnalyses directory. Use `]` to enter Pkg mode. 
```julia

(@v1.11) pkg> activate .
  Activating project at `~/Projects/EntrapmentAnalysesJmod/EntrapmentAnalyses`

julia> using Revise, EntrapmentAnalyses
Precompiling EntrapmentAnalyses...
  1 dependency successfully precompiled in 6 seconds. 273 already precompiled.
```


## Input Data Requirements

### PSM Results (Parquet Format)

Your PSM results file should contain the following columns:

- `stripped_seq`: Peptide sequence without modifications (String)
- `z`: Precursor charge state (Integer)
- `PredVal`: Prediction score from your search engine (Float)
- `decoy`: Boolean flag indicating decoy status (Bool)
- `file_name`: Identifier for the source file (String)
- `protein`: Protein identifiers, semicolon-separated for multi-mapping peptides (String)

### Spectral Library (TSV Format)

The spectral library should contain:

- `PeptideSequence`: Modified peptide sequence (String)
- `PrecursorCharge`: Charge state (Integer)
- `EntrapmentGroupId`: 0 for original peptides, >0 for entrapments (Integer)
- `PrecursorIdx`: Unique identifier linking original/entrapment pairs (Integer)

## Basic Usage

### Precursor-Level Analysis

```julia
using EntrapmentAnalyses

# Single file analysis
results = run_efdr_analysis(
    ["data/search_results.parquet"],
    "data/spectral_library.tsv"
)

# Multiple file analysis
results = run_efdr_analysis(
    ["data/run1.parquet", "data/run2.parquet", "data/run3.parquet"],
    "data/spectral_library.tsv";
    output_dir="results/efdr_analysis"
)
```

### Protein-Level Analysis

```julia
# Protein-level EFDR analysis
protein_results = run_protein_efdr_analysis(
    ["data/search_results.parquet"],
    "data/spectral_library.tsv";
    output_dir="results/protein_analysis",
    protein_q_threshold=0.01  # 1% protein-level FDR
)
```

## Advanced Options

### Customizing EFDR Methods

```julia
# Run only combined EFDR
results = run_efdr_analysis(
    parquet_files,
    library_file;
    efdr_method=:combined
)

# Run only paired EFDR
results = run_efdr_analysis(
    parquet_files,
    library_file;
    efdr_method=:paired
)

# Run both methods (default)
results = run_efdr_analysis(
    parquet_files,
    library_file;
    efdr_method=:both
)
```

### Custom Output Options

```julia
results = run_efdr_analysis(
    parquet_files,
    library_file;
    output_dir="custom_output",
    save_plots=true,           # Generate and save plots
    save_tables=true,          # Save result tables
    generate_report=true       # Create markdown report
)
```

## Working with Results

### Understanding the Output

The analysis generates several output files:

1. **EFDR Results Table** (`efdr_results.tsv`): Contains EFDR values at different q-value thresholds
2. **Plots** (`plots/` directory): EFDR vs q-value plots for each method
3. **Analysis Report** (`analysis_report.md`): Comprehensive markdown report with embedded plots

### Accessing Results Programmatically

```julia
# The function returns a DataFrame with results
results = run_efdr_analysis(files, library)

# Access EFDR values
combined_efdr = results[results.method .== "combined", :]
paired_efdr = results[results.method .== "paired", :]

# Plot custom visualizations
using Plots
plot(combined_efdr.q_value, combined_efdr.efdr, 
     label="Combined EFDR", 
     xlabel="Q-value", 
     ylabel="EFDR")
```

## Troubleshooting

### Common Issues

1. **Missing Values in Data**
   - The package automatically handles missing values
   - Check console output for warnings about replaced values

2. **Memory Issues with Large Datasets**
   - Process files in smaller batches
   - Reduce `local_qval_threshold` to filter more aggressively

3. **Pairing Warnings**
   - "No complement found" warnings are normal for unpaired peptides
   - Ensure your spectral library has correct `PrecursorIdx` values

### Performance Tips

1. **Use Local SSD Storage**: Reading from fast local storage improves performance
2. **Batch Processing**: For many files, process in batches of 10-20
3. **Pre-filter Data**: Remove low-quality PSMs before analysis

## Best Practices

1. **Q-value Thresholds**: Start with default 5% threshold, adjust based on your needs
2. **Multiple Runs**: Always analyze multiple technical replicates together
3. **Validation**: Compare EFDR results with traditional target-decoy FDR
4. **Documentation**: Keep track of analysis parameters for reproducibility