# Tutorial

This tutorial walks through a complete analysis workflow using EntrapmentAnalyses.jl.

## Setup

First, let's set up our environment and load the necessary packages:

```julia
using EntrapmentAnalyses
using DataFrames
using CSV
```

## Example 1: Basic Precursor-Level Analysis

Let's start with a simple analysis of a single proteomics run.

### Step 1: Prepare Your Data

Ensure you have:
1. A PSM results file in Parquet format
2. A spectral library in TSV format

```julia
# Define file paths
psm_file = "data/example_run.parquet"
library_file = "data/spectral_library.tsv"
output_dir = "tutorial_output"
```

### Step 2: Run the Analysis

```julia
# Run EFDR analysis with default parameters
results = run_efdr_analysis(
    [psm_file],
    library_file;
    output_dir=output_dir
)

# View the results
println("EFDR Results:")
println(results)
```

### Step 3: Interpret the Results

The results DataFrame contains:
- `q_value`: The q-value threshold
- `efdr`: The empirical FDR at that threshold
- `method`: Either "combined" or "paired"
- `run`: The file name (for multi-file analyses)

## Example 2: Multi-Run Analysis

When analyzing multiple runs, the package automatically handles per-file processing:

```julia
# Define multiple PSM files
psm_files = [
    "data/run1.parquet",
    "data/run2.parquet",
    "data/run3.parquet"
]

# Run analysis across all files
multi_results = run_efdr_analysis(
    psm_files,
    library_file;
    output_dir="multi_run_output",
    efdr_method=:both  # Calculate both combined and paired EFDR
)

# View results for each run
for run in unique(multi_results.run)
    run_data = filter(row -> row.run == run, multi_results)
    println("\nResults for $run:")
    println(run_data)
end
```

## Example 3: Protein-Level Analysis

For protein-level FDR analysis:

```julia
# Run protein-level analysis
protein_results = run_protein_efdr_analysis(
    psm_files,
    library_file;
    output_dir="protein_output",
    protein_q_threshold=0.01,  # 1% protein FDR
    efdr_method=:paired       # Use paired EFDR for proteins
)

# The function performs protein rollup automatically
println("Protein-level EFDR results saved to: protein_output/")
```

## Example 4: Custom Analysis Pipeline

For more control, you can use the lower-level API:

```julia
# Load data manually
psm_df = load_parquet_results(psm_file)
spec_lib_df = load_spectral_library(library_file)

# Add pairing information
compute_pairing_vectors!(psm_df, spec_lib_df)

# Calculate q-values per file
calculate_qvalues_per_file!(psm_df, :PredVal)

# Filter by q-value
filtered_df = filter(row -> row.q_value <= 0.05, psm_df)

# Calculate EFDR
combined_efdr_df = calculate_combined_efdr(filtered_df)
paired_efdr_df = calculate_paired_efdr(filtered_df)

# Create custom plots
using Plots
plot(combined_efdr_df.q_value, combined_efdr_df.efdr,
     label="Combined EFDR",
     xlabel="Q-value",
     ylabel="EFDR",
     linewidth=2)
```

## Example 5: Handling Special Cases

### Working with Missing Values

The package automatically handles missing values, but you can check for them:

```julia
# Load data
psm_df = load_parquet_results(psm_file)

# Check for missing values
for col in names(psm_df)
    n_missing = sum(ismissing.(psm_df[!, col]))
    if n_missing > 0
        println("Column $col has $n_missing missing values")
    end
end
```

### Debugging Pairing Issues

If you encounter pairing warnings:

```julia
# Get detailed pairing information
psm_df = load_parquet_results(psm_file)
spec_lib_df = load_spectral_library(library_file)

# Add pairing with verbose output
compute_pairing_vectors!(psm_df, spec_lib_df)

# Check pairing statistics
n_paired = sum(psm_df.pair_id .!= -1)
n_unpaired = sum(psm_df.pair_id .== -1)
println("Paired peptides: $n_paired")
println("Unpaired peptides: $n_unpaired")
```

## Example 6: Batch Processing

For large-scale analyses:

```julia
# Process files in batches
all_files = readdir("data/psm_results", join=true)
batch_size = 10

for i in 1:batch_size:length(all_files)
    batch_end = min(i + batch_size - 1, length(all_files))
    batch_files = all_files[i:batch_end]
    
    println("Processing batch $(div(i, batch_size) + 1)")
    
    results = run_efdr_analysis(
        batch_files,
        library_file;
        output_dir="batch_$(div(i, batch_size) + 1)"
    )
end
```

## Tips and Best Practices

1. **Start Small**: Test your analysis on a single file before processing large batches
2. **Monitor Progress**: The package shows progress bars for long operations
3. **Check Outputs**: Always review the generated plots and reports
4. **Save Intermediate Results**: For complex analyses, save intermediate DataFrames

## Next Steps

- Explore the [API Reference](@ref) for detailed function documentation
- Read about [advanced features](guide.md) in the User Guide
- Check the [GitHub repository](https://github.com/nathanwamsley/EntrapmentAnalyses.jl) for updates