# API Reference

```@meta
CurrentModule = EntrapmentAnalyses
```

## Main API Functions

This is the primary function for running EFDR analyses:

```@docs
run_efdr_analysis
```

## Data Loading

Functions for loading input data:

```@docs
load_parquet
load_spectral_library
```

## EFDR Calculation

Core EFDR calculation method:

```@docs
calculate_efdr
```

## Q-value Calculation

Functions for q-value computation:

```@docs
calculate_qvalues!
calculate_qvalues_per_file!
calculate_global_qvalues!
monotonize!
```

## Pairing System

Functions for computing peptide pairings:

```@docs
compute_pairing_vectors!
```

## Visualization

Functions for generating plots and reports:

```@docs
plot_efdr_comparison
plot_combined_efdr
plot_paired_efdr
plot_efdr_comparison_both_methods
generate_analysis_report
```

## Type Definitions

Key types used throughout the package:

```@docs
EFDRMethod
CombinedEFDR
PairedEFDR
```

## Index

```@index
```