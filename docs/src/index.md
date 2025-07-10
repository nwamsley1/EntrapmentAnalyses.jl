```@meta
CurrentModule = EntrapmentAnalyses
```

# EntrapmentAnalyses.jl

[![Build Status](https://github.com/nwamsley1/EntrapmentAnalyses.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/nwamsley1/EntrapmentAnalyses.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![codecov](https://codecov.io/gh/nwamsley1/EntrapmentAnalyses.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/nwamsley1/EntrapmentAnalyses.jl)

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

## Documentation Contents

```@contents
Pages = ["guide.md", "tutorial.md", "api.md"]
Depth = 2
```

### Workflow

1. Load PSM results and spectral library data
2. Compute peptide pairings with plex-aware complement scoring
3. Calculate per-file q-values
4. Perform EFDR analysis (combined and/or paired)
5. Generate visualizations and reports
6. Optional: Perform protein-level analysis
