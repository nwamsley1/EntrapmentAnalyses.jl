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

# Export data loading functions
export load_parquet, load_spectral_library

# Export main API functions
export run_efdr_analysis, run_protein_efdr_analysis
export analyze_combined_efdr, analyze_paired_efdr, analyze_proteins

# Export types
export EFDRMethod, CombinedEFDR, PairedEFDR

# Export utility functions
export calculate_efdr, monotonize!
export compute_pairing_vectors, compute_pairing_vectors!
export calculate_qvalues!, calculate_global_qvalues!, calculate_qvalues_per_file!

# Export plotting functions
export plot_efdr_comparison, plot_protein_comparison
export plot_combined_efdr, plot_paired_efdr, plot_efdr_comparison_both_methods
export generate_analysis_report

end # module EntrapmentAnalyses
