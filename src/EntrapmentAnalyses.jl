module EntrapmentAnalyses

using DataFrames
using CSV
using Parquet2: Dataset
using Plots
using Printf
using Statistics
using Dates
using ProgressBars

# Core components - order matters for dependencies
include("io/data_loaders.jl")
include("core/pairing.jl")
include("core/efdr_methods.jl")
include("core/scoring.jl")
include("analysis/qvalue_calculation.jl")
include("analysis/protein_analysis.jl")
include("analysis/efdr_analysis.jl")
include("plotting/visualization.jl")
include("api.jl")

# Export main API functions
export run_efdr_analysis, run_protein_efdr_analysis

# Export types
export EFDRMethod, CombinedEFDR, PairedEFDR

# Export utility functions
export calculate_efdr, monotonize!
export compute_pairing_vectors
export calculate_qvalues!, calculate_global_qvalues!

end # module EntrapmentAnalyses
