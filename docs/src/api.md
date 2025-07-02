# API Reference

```@meta
CurrentModule = EntrapmentAnalyses
```

## Main API Functions

These are the primary functions for running EFDR analyses:

```@docs
run_efdr_analysis
run_protein_efdr_analysis
```

## Data Loading

Functions for loading input data:

```@docs
load_parquet_results
load_spectral_library
```

## Pairing System

Functions for computing peptide pairings:

```@docs
compute_pairing_vectors!
init_entrapment_pairs_dict
add_plex_complement_scores!
```

## EFDR Calculation

Core EFDR calculation methods:

```@docs
calculate_combined_efdr
calculate_paired_efdr
```

## Q-value Calculation

Functions for q-value computation:

```@docs
calculate_qvalues_per_file!
calculate_qvalues!
monotonize!
```

## Protein Analysis

Functions for protein-level analysis:

```@docs
perform_protein_rollup
calculate_protein_efdr_per_run
```

## Visualization

Functions for generating plots and reports:

```@docs
plot_efdr_results
plot_efdr_comparison
generate_analysis_report
```

## Utility Functions

Helper functions and utilities:

```@docs
handle_missing_values!
strip_modifications
```

## Type Definitions

Key types used throughout the package:

```@docs
PeptideKey
PlexPairKey
ScorePair
```

## Index

```@index
```