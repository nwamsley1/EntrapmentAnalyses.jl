# Implementation Plan: EntrapmentAnalyses.jl Module (Revised)

## Executive Summary
This plan outlines the development of the EntrapmentAnalyses.jl module that processes proteomics data in Parquet/TSV format with an efficient, typed implementation of empirical FDR calculation using entrapment sequences.

## Core Design Principles
1. **Native Format**: Accept Parquet (PSM results) and TSV (spectral library) directly
2. **Performance**: Pre-compute pairing information to avoid dictionary lookups in main loops
3. **Type Safety**: Use typed structs (CombinedEFDR, PairedEFDR) for EFDR methods
4. **Clean API**: No references to implementation origins in naming

## Phase 1: Core Module Structure

### 1.1 Module Definition
```julia
# src/EntrapmentAnalyses.jl
module EntrapmentAnalyses

using DataFrames
using CSV
using Parquet2
using Plots
using Printf
using Statistics
using Dates

# Core components
include("io/data_loaders.jl")
include("core/pairing.jl")
include("core/efdr_methods.jl")
include("core/scoring.jl")
include("analysis/efdr_analysis.jl")
include("analysis/qvalue_calculation.jl")
include("plotting/visualization.jl")
include("api.jl")

# Export main API
export run_efdr_analysis

# Export types
export EFDRMethod, CombinedEFDR, PairedEFDR

# Export core functions
export calculate_efdr, monotonize!
export compute_pairing_vectors

end
```

## Phase 2: Efficient Pairing System

### 2.1 Pre-computed Pairing Vectors
```julia
# src/core/pairing.jl

"""
Pre-compute pairing information as vectors for efficient access
Returns vectors that can be used without dictionary lookups
"""
function compute_pairing_vectors(
    library_df::DataFrame,
    results_df::DataFrame;
    seq_col::Symbol = :PeptideSequence,
    charge_col::Symbol = :PrecursorCharge,
    entrap_col::Symbol = :EntrapmentGroupId,
    pair_col::Symbol = :PrecursorIdx
)
    # Create lookup vectors indexed by result row
    n_results = nrow(results_df)
    
    # Pre-allocate vectors
    is_original = Vector{Bool}(undef, n_results)
    pair_indices = Vector{Int}(undef, n_results)
    entrap_labels = Vector{Int}(undef, n_results)
    
    # Build temporary lookup for O(1) access during setup
    temp_lookup = Dict{Tuple{String, Int}, Int}()
    for i in 1:nrow(library_df)
        key = (library_df[i, seq_col], library_df[i, charge_col])
        temp_lookup[key] = i
    end
    
    # Fill vectors
    for i in 1:n_results
        seq = results_df[i, :stripped_seq]
        z = results_df[i, :z]
        lib_idx = temp_lookup[(seq, z)]
        
        is_original[i] = library_df[lib_idx, entrap_col] == 0
        pair_indices[i] = library_df[lib_idx, pair_col]
        entrap_labels[i] = library_df[lib_idx, entrap_col]
    end
    
    # Find complement indices for each row
    complement_indices = Vector{Int}(undef, n_results)
    
    # Group by pair_id for efficient lookup
    for pair_id in unique(pair_indices)
        mask = pair_indices .== pair_id
        indices = findall(mask)
        
        for idx in indices
            if is_original[idx]
                # Find entrapment partner
                for other_idx in indices
                    if !is_original[other_idx]
                        complement_indices[idx] = other_idx
                        break
                    end
                end
            else
                # Find original partner
                for other_idx in indices
                    if is_original[other_idx]
                        complement_indices[idx] = other_idx
                        break
                    end
                end
            end
        end
    end
    
    return (
        is_original = is_original,
        pair_indices = pair_indices,
        entrap_labels = entrap_labels,
        complement_indices = complement_indices
    )
end
```

## Phase 3: EFDR Methods with Pre-computed Data

### 3.1 Typed EFDR Structs
```julia
# src/core/efdr_methods.jl

abstract type EFDRMethod end

struct CombinedEFDR{T<:Real} <: EFDRMethod
    scores::Vector{T}
    entrap_labels::Vector{Int}
    qvals::Vector{T}
    r::T
end

struct PairedEFDR{T<:Real} <: EFDRMethod
    scores::Vector{T}
    complement_scores::Vector{T}
    is_original::Vector{Bool}
    entrap_labels::Vector{Int}
    qvals::Vector{T}
    pair_indices::Vector{Int}
    r::T
end

"""
Calculate empirical FDR using pre-computed vectors
"""
function calculate_efdr(method::PairedEFDR)
    n = length(method.scores)
    efdr = zeros(eltype(method.qvals), n)
    
    # Sort by q-value (ascending), then score (descending)
    sort_order = sortperm(collect(zip(method.qvals, -method.scores)))
    
    # Pre-compute score relationships for all pairs
    # This avoids repeated lookups in the main loop
    score_relationships = Vector{Symbol}(undef, n)
    
    for i in 1:n
        if method.is_original[i]
            score_relationships[i] = :original
        else
            orig_score = method.complement_scores[i]
            entr_score = method.scores[i]
            
            if entr_score > orig_score
                score_relationships[i] = :entrap_wins
            elseif entr_score == orig_score
                score_relationships[i] = :tie
            else
                score_relationships[i] = :original_wins
            end
        end
    end
    
    # Main calculation loop - no dictionary lookups
    for i in 1:n
        Nτ, Nε, Nεsτ, Nετs = 0, 0, 0, 0
        s = method.scores[sort_order[i]]
        
        for j in 1:i
            idx = sort_order[j]
            
            if method.is_original[idx]
                Nτ += 1
            else
                Nε += 1
                
                # Use pre-computed relationships
                if score_relationships[idx] == :entrap_wins
                    orig_score = method.complement_scores[idx]
                    if orig_score >= s
                        Nετs += 1
                    end
                elseif method.scores[idx] >= s
                    orig_score = method.complement_scores[idx]
                    if orig_score < s
                        Nεsτ += 1
                    end
                end
            end
        end
        
        if Nε + Nτ > 0
            efdr[sort_order[i]] = min(1.0, (Nε + Nεsτ + 2*Nετs) / (Nε + Nτ))
        end
    end
    
    return efdr
end
```

### 3.2 Helper Functions
```julia
"""
Monotonize FDR values to ensure non-decreasing property
"""
function monotonize!(values::AbstractVector{T}) where T<:AbstractFloat
    current_min = T(1.0)
    
    # Work backwards through sorted results
    for i in length(values):-1:1
        if values[i] > current_min
            values[i] = current_min
        else
            current_min = values[i]
        end
    end
    
    return values
end
```

## Phase 4: Analysis Pipeline

### 4.1 Main API Function
```julia
# src/api.jl

"""
Run empirical FDR analysis on proteomics data
"""
function run_efdr_analysis(
    parquet_files::Vector{String},
    library_path::String;
    output_dir::String = "efdr_output",
    score_col::Symbol = :PredVal,
    methods::Vector{DataType} = [CombinedEFDR, PairedEFDR],
    global_qval_threshold::Float64 = 0.01,
    r_lib::Float64 = 1.0,
    analyze_proteins::Bool = true
)
    # Step 1: Load data
    println("Loading PSM results...")
    results_df = load_parquet_results(parquet_files)
    
    println("Loading spectral library...")
    library_df = DataFrame(CSV.File(library_path))
    
    # Step 2: Add entrapment information
    println("Computing pairing vectors...")
    pairing_info = compute_pairing_vectors(library_df, results_df)
    
    # Step 3: Filter decoys and calculate q-values
    println("Calculating q-values...")
    filter!(:decoy => ==(false), results_df)
    calculate_qvalues!(results_df, score_col)
    
    # Step 4: Prepare scores for EFDR calculation
    scores = Float64.(results_df[!, score_col])
    qvals = Float64.(results_df[!, :local_qvalue])
    
    # Get complement scores efficiently
    complement_scores = Vector{Float64}(undef, nrow(results_df))
    for i in 1:nrow(results_df)
        if pairing_info.complement_indices[i] > 0
            complement_scores[i] = scores[pairing_info.complement_indices[i]]
        else
            complement_scores[i] = 0.0  # No pair found
        end
    end
    
    # Step 5: Calculate EFDR for each method
    efdr_results = Dict{Symbol, Vector{Float64}}()
    
    for method_type in methods
        println("Calculating $(method_type) EFDR...")
        
        if method_type == CombinedEFDR
            method = CombinedEFDR(
                scores,
                pairing_info.entrap_labels,
                qvals,
                r_lib
            )
            efdr_results[:combined_efdr] = calculate_efdr(method)
            
        elseif method_type == PairedEFDR
            method = PairedEFDR(
                scores,
                complement_scores,
                pairing_info.is_original,
                pairing_info.entrap_labels,
                qvals,
                pairing_info.pair_indices,
                r_lib
            )
            efdr_results[:paired_efdr] = calculate_efdr(method)
        end
        
        # Monotonize results
        monotonize!(efdr_results[Symbol(lowercase(string(method_type)))])
    end
    
    # Step 6: Add EFDR columns to results
    for (col_name, efdr_vals) in efdr_results
        results_df[!, col_name] = efdr_vals
    end
    
    # Step 7: Generate outputs
    mkpath(output_dir)
    
    # Save results
    CSV.write(joinpath(output_dir, "results_with_efdr.tsv"), results_df, delim='\t')
    
    # Generate plots
    for (col_name, efdr_vals) in efdr_results
        plot_efdr_comparison(
            results_df[!, :local_qvalue],
            efdr_vals,
            title = "EFDR Analysis - $(col_name)",
            output_path = joinpath(output_dir, "efdr_comparison_$(col_name).pdf")
        )
    end
    
    return results_df
end
```

## Phase 5: Data Loading

### 5.1 Parquet Loader
```julia
# src/io/data_loaders.jl

function load_parquet_results(filepaths::Vector{String})
    dfs = DataFrame[]
    
    for filepath in filepaths
        push!(dfs, DataFrame(Parquet2.Dataset(filepath); copycols=true))
    end
    
    results_df = vcat(dfs...)
    
    # Add entrapment group based on sequence matching
    # (This would be done after loading library)
    
    return results_df
end
```

## Phase 6: Visualization

### 6.1 Plotting Functions
```julia
# src/plotting/visualization.jl

function plot_efdr_comparison(
    qvals::Vector{Float64},
    efdr_vals::Vector{Float64};
    title::String = "EFDR Comparison",
    output_path::String = "efdr_comparison.pdf"
)
    p = plot(
        size=(600, 450),
        xlim = (0, 0.05), 
        ylim = (0, 0.05),
        xlabel = "FDR", 
        ylabel = "Entrapment FDR",
        title = title,
        titlefontsize = 16,
        guidefontsize = 16,
        tickfontsize = 12
    )
    
    # Diagonal reference line
    plot!(p, [0, 0.05], [0, 0.05], 
          color = :black, lw = 3, ls = :dash, label = nothing)
    
    # EFDR comparison
    plot!(p, qvals, efdr_vals,
          lw = 3,
          color = RGB(0.39215686, 0.58431373, 0.92941176),
          alpha = 0.75,
          label = nothing)
    
    savefig(p, output_path)
    
    return p
end
```

## Key Improvements Over Original Plan

1. **No Dictionary Lookups**: Pre-compute all pairing information as vectors
2. **Clean Naming**: No references to "notebook" or "module"
3. **Direct Format Support**: Native Parquet/TSV handling without conversion
4. **Efficient Algorithm**: O(n²) but with minimal overhead per iteration
5. **Type Safety**: Fully typed structs for EFDR methods

## Implementation Timeline
- **Day 1**: Core structure and data loading
- **Day 2**: Efficient pairing system
- **Day 3**: EFDR methods implementation
- **Day 4**: Analysis pipeline and API
- **Day 5**: Visualization and testing

## Performance Considerations
- Pre-computed vectors eliminate dictionary overhead
- Type-stable operations throughout
- Minimal allocations in inner loops
- Efficient sorting with sortperm