# User Guide

This guide provides comprehensive instructions for using EntrapmentAnalyses.jl in your proteomics workflows.

## Installation

EntrapmentAnalyses.jl requires Julia 1.10 or later. Navigate to the EntrapmentAnalyses directory. Use `]` to enter Pkg mode. 

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
- `channel`: Integer unique to a tag plex 
- `file_name`: Identifier for the source file (String)

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
