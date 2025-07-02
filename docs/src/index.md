```@meta
CurrentModule = EntrapmentAnalyses
```

# EntrapmentAnalyses.jl

[![CI](https://github.com/nathanwamsley/EntrapmentAnalyses.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/nathanwamsley/EntrapmentAnalyses.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/nathanwamsley/EntrapmentAnalyses.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/nathanwamsley/EntrapmentAnalyses.jl)

EntrapmentAnalyses.jl is a Julia package for entrapment-based false discovery rate (FDR) analysis in proteomics data. It implements both combined and paired empirical FDR (EFDR) methods for precursor-level and protein-level analyses.

## Features

- **Data Loading**: Support for Parquet and TSV file formats
- **EFDR Calculation**: Both combined and paired empirical FDR methods
- **Protein Analysis**: Protein-level rollup and per-run analysis
- **Visualization**: Automated generation of EFDR plots and reports
- **Performance**: Optimized for large-scale proteomics datasets

## Installation

```julia
using Pkg
Pkg.add("EntrapmentAnalyses")
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

# Run protein-level EFDR analysis
protein_results = run_protein_efdr_analysis(
    ["data/psm_results.parquet"],
    "data/spectral_library.tsv";
    output_dir="output"
)
```

## Documentation Contents

```@contents
Pages = ["guide.md", "tutorial.md", "api.md"]
Depth = 2
```

## Package Overview

EntrapmentAnalyses.jl provides a comprehensive toolkit for analyzing entrapment-based FDR in proteomics experiments. The package is designed to handle the complete workflow from data loading through analysis to visualization.

### Key Components

1. **Data Loading**: Robust handling of various input formats with automatic missing value management
2. **Pairing System**: Sophisticated peptide pairing that respects both file and plex boundaries
3. **EFDR Methods**: Implementation of both combined and paired EFDR calculations
4. **Protein Analysis**: Tools for protein-level rollup and per-run analysis
5. **Visualization**: Automated plot generation with customizable parameters

### Workflow

1. Load PSM results and spectral library data
2. Compute peptide pairings with plex-aware complement scoring
3. Calculate per-file q-values
4. Perform EFDR analysis (combined and/or paired)
5. Generate visualizations and reports
6. Optional: Perform protein-level analysis

## Getting Help

- Check the [User Guide](@ref) for detailed usage instructions
- Follow the [Tutorial](@ref) for a step-by-step walkthrough
- Consult the [API Reference](@ref) for function documentation
- Report issues on [GitHub](https://github.com/nathanwamsley/EntrapmentAnalyses.jl/issues)