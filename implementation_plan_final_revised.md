# Implementation Plan: EntrapmentAnalyses.jl Module (Final Revised Version)

## Reference Implementations

### Jupyter Notebook Reference
- **Path**: `/Users/nathanwamsley/Data/May2025/kmd_jmod_search/9plex/05202025/entrapment_test.ipynb`
- **Description**: Original implementation using dictionary-based pairing for multiplexed proteomics data
- **Key Functions**: 
  - `getPairedEmpiricalFdr()` - Main EFDR calculation
  - `initEntrapPairsDict()` - Creates pairing dictionaries
  - `monotonize!()` - Ensures FDR monotonicity
- **Data Format**: Parquet files for PSM results, TSV for spectral library

### Pioneer.jl Module Reference
- **Path**: `/Users/nathanwamsley/Projects/Pioneer.jl/test/entrapment_analyses/`
- **Description**: Structured module implementation with typed EFDR methods
- **Key Features**:
  - Abstract type system for EFDR methods
  - Arrow format I/O
  - Comprehensive reporting

## Module Structure

```
/Users/nathanwamsley/Projects/EntrapmentAnalysesJmod/EntrapmentAnalyses/
├── Project.toml
├── src/
│   ├── EntrapmentAnalyses.jl       # Main module file
│   ├── io/
│   │   └── data_loaders.jl         # Parquet/TSV loading
│   ├── core/
│   │   ├── pairing.jl              # Efficient pairing vectors
│   │   ├── efdr_methods.jl         # EFDR calculation implementations
│   │   └── scoring.jl              # Score utilities
│   ├── analysis/
│   │   ├── efdr_analysis.jl        # Analysis coordination
│   │   ├── qvalue_calculation.jl   # Q-value calculations
│   │   └── protein_analysis.jl     # Protein-level rollup
│   ├── plotting/
│   │   └── visualization.jl         # Plot generation
│   └── api.jl                      # Public API functions
└── test/
    ├── runtests.jl
    ├── unit/
    │   ├── test_pairing.jl
    │   ├── test_efdr_methods.jl
    │   └── test_qvalue.jl
    └── integration/
        ├── test_full_pipeline.jl
        └── test_data/
            ├── sample_results.parquet
            └── sample_library.tsv
```

## Phase 1: Core Module Setup with Dependencies

### 1.1 Project.toml Configuration
```toml
name = "EntrapmentAnalyses"
uuid = "60062a3e-e55a-4920-8913-bf16cdadf465"
authors = ["Nathan Wamsley <wamsleynathan@gmail.com>"]
version = "0.1.0"

[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
Dates = "ade2ca70-3891-5945-98fb-dc099432e06a"
Parquet2 = "98572fba-bba0-415d-956f-fa77e587d26d"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"
ProgressBars = "49802e3a-d2f1-5c88-81d8-b72133a6f568"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[compat]
CSV = "0.10"
DataFrames = "1.6"
Parquet2 = "0.2"
Plots = "1.40"
ProgressBars = "1.5"
julia = "1.9"
```

### 1.2 Main Module File
```julia
# src/EntrapmentAnalyses.jl
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

end # module
```

## Phase 2: Data Loading with Validation

### 2.1 Data Loaders Implementation
```julia
# src/io/data_loaders.jl

"""
    load_parquet_results(filepaths::Vector{String})

Load and combine multiple Parquet files containing PSM results.

# Expected columns:
- stripped_seq: String - Peptide sequence without modifications
- z: Integer - Charge state
- PredVal: Float - Prediction score (higher is better)
- decoy: Boolean - Decoy status
- channel: Integer - Multiplex channel (optional, will add dummy if missing)
- file_name: String - Source file identifier
- protein: String - Protein group (for protein-level analysis)

# Returns
- DataFrame with combined results
"""
function load_parquet_results(filepaths::Vector{String})
    # Validate files exist
    for filepath in filepaths
        if !isfile(filepath)
            error("File not found: $filepath")
        end
    end
    
    # Load all files
    dfs = DataFrame[]
    
    for (i, filepath) in enumerate(filepaths)
        println("Loading file $i/$(length(filepaths)): $(basename(filepath))")
        df = DataFrame(Dataset(filepath); copycols=true)
        push!(dfs, df)
    end
    
    # Combine dataframes
    results_df = vcat(dfs...)
    
    # Validate required columns
    required_cols = [:stripped_seq, :z, :PredVal, :decoy]
    missing_cols = setdiff(required_cols, names(results_df, Symbol))
    
    if !isempty(missing_cols)
        error("Missing required columns: $(join(missing_cols, ", "))")
    end
    
    # Add dummy channel if missing
    if !hasproperty(results_df, :channel)
        println("No channel column detected, adding dummy channel 0")
        results_df[!, :channel] = zero(UInt8)
    end
    
    println("Loaded $(nrow(results_df)) total PSMs")
    
    return results_df
end

"""
    load_spectral_library(filepath::String)

Load spectral library from TSV file.

# Expected columns:
- PeptideSequence: String - Modified peptide sequence
- PrecursorCharge: Integer - Charge state
- EntrapmentGroupId: Integer - 0 for originals, >0 for entrapments
- PrecursorIdx: Integer - Unique pair identifier

# Returns
- DataFrame with library information
"""
function load_spectral_library(filepath::String)
    if !isfile(filepath)
        error("Library file not found: $filepath")
    end
    
    println("Loading spectral library from: $(basename(filepath))")
    library_df = DataFrame(CSV.File(filepath))
    
    # Validate required columns
    required_cols = [:PeptideSequence, :PrecursorCharge, :EntrapmentGroupId, :PrecursorIdx]
    missing_cols = setdiff(required_cols, names(library_df, Symbol))
    
    if !isempty(missing_cols)
        error("Missing required columns in library: $(join(missing_cols, ", "))")
    end
    
    # Summary statistics
    n_targets = sum(library_df.EntrapmentGroupId .== 0)
    n_entrapments = sum(library_df.EntrapmentGroupId .> 0)
    
    println("Library loaded: $n_targets targets, $n_entrapments entrapments")
    
    return library_df
end
```

## Phase 3: Efficient Pairing System with Progress Monitoring

### 3.1 Pairing Vector Computation
```julia
# src/core/pairing.jl

"""
    compute_pairing_vectors(library_df::DataFrame, results_df::DataFrame; kwargs...)

Pre-compute pairing information as vectors for efficient EFDR calculation.
Avoids dictionary lookups in main computation loops.

# Arguments
- `library_df`: Spectral library with entrapment information
- `results_df`: PSM results to analyze

# Keyword Arguments
- `seq_col`: Column name for sequences in library (default: :PeptideSequence)
- `charge_col`: Column name for charge in library (default: :PrecursorCharge)
- `entrap_col`: Column name for entrapment group (default: :EntrapmentGroupId)
- `pair_col`: Column name for pair ID (default: :PrecursorIdx)
- `show_progress`: Show progress bar (default: true)

# Returns
Named tuple with:
- `is_original`: Bool vector indicating if each result is an original sequence
- `pair_indices`: Int vector with pair IDs for each result
- `entrap_labels`: Int vector with entrapment group IDs
- `complement_indices`: Int vector with indices of paired sequences
"""
function compute_pairing_vectors(
    library_df::DataFrame,
    results_df::DataFrame;
    seq_col::Symbol = :PeptideSequence,
    charge_col::Symbol = :PrecursorCharge,
    entrap_col::Symbol = :EntrapmentGroupId,
    pair_col::Symbol = :PrecursorIdx,
    show_progress::Bool = true
)
    n_results = nrow(results_df)
    println("Computing pairing vectors for $n_results PSMs...")
    
    # Pre-allocate output vectors
    is_original = Vector{Bool}(undef, n_results)
    pair_indices = Vector{Int}(undef, n_results)
    entrap_labels = Vector{Int}(undef, n_results)
    complement_indices = fill(-1, n_results)  # -1 indicates no pair found
    
    # Build efficient lookup structure
    # Key: (sequence, charge) -> Value: library row index
    lib_lookup = Dict{Tuple{String, Int}, Int}()
    
    for i in 1:nrow(library_df)
        key = (library_df[i, seq_col], library_df[i, charge_col])
        lib_lookup[key] = i
    end
    
    # Fill basic pairing information
    pb = show_progress ? ProgressBar(1:n_results) : 1:n_results
    set_description(pb, "Mapping sequences to library")
    
    for i in pb
        seq = results_df[i, :stripped_seq]
        z = results_df[i, :z]
        key = (seq, z)
        
        if haskey(lib_lookup, key)
            lib_idx = lib_lookup[key]
            is_original[i] = library_df[lib_idx, entrap_col] == 0
            pair_indices[i] = library_df[lib_idx, pair_col]
            entrap_labels[i] = library_df[lib_idx, entrap_col]
        else
            error("Sequence not found in library: $seq with charge $z")
        end
    end
    
    # Build pair lookup for complement finding
    # pair_id -> [result_indices...]
    pair_groups = Dict{Int, Vector{Int}}()
    
    for i in 1:n_results
        pid = pair_indices[i]
        if !haskey(pair_groups, pid)
            pair_groups[pid] = Int[]
        end
        push!(pair_groups[pid], i)
    end
    
    # Find complement indices
    println("Finding complement sequences...")
    pb2 = show_progress ? ProgressBar(1:n_results) : 1:n_results
    set_description(pb2, "Pairing complements")
    
    for i in pb2
        pid = pair_indices[i]
        group = pair_groups[pid]
        
        # Find complement within the same pair group
        for j in group
            if i != j && is_original[i] != is_original[j]
                complement_indices[i] = j
                break
            end
        end
    end
    
    # Validate pairing
    n_unpaired = sum(complement_indices .== -1)
    if n_unpaired > 0
        @warn "$n_unpaired sequences have no complement pair"
    end
    
    return (
        is_original = is_original,
        pair_indices = pair_indices,
        entrap_labels = entrap_labels,
        complement_indices = complement_indices
    )
end
```

## Phase 4: EFDR Methods with Progress Monitoring

### 4.1 EFDR Method Implementations
```julia
# src/core/efdr_methods.jl

abstract type EFDRMethod end

"""
Combined EFDR method - standard empirical FDR calculation
"""
struct CombinedEFDR{T<:Real} <: EFDRMethod
    scores::Vector{T}
    entrap_labels::Vector{Int}
    qvals::Vector{T}
    r::T
end

"""
Paired EFDR method - considers score relationships between pairs
"""
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
    calculate_efdr(method::CombinedEFDR; show_progress=true)

Calculate combined empirical FDR.
"""
function calculate_efdr(method::CombinedEFDR; show_progress::Bool=true)
    n = length(method.scores)
    efdr = zeros(eltype(method.qvals), n)
    
    # Validate sort order
    if !issorted(method.qvals)
        @warn "Q-values are not sorted. Results may be incorrect."
    end
    
    # Sort by q-value (ascending), then score (descending)
    sort_order = sortperm(collect(zip(method.qvals, -method.scores)))
    
    Nτ, Nε = 0, 0
    
    pb = show_progress ? ProgressBar(1:n) : 1:n
    set_description(pb, "Calculating Combined EFDR")
    
    for i in pb
        idx = sort_order[i]
        
        if method.entrap_labels[idx] == 0
            Nτ += 1
        else
            Nε += 1
        end
        
        if Nτ + Nε > 0
            efdr[idx] = min(1.0, (Nε * (1 + 1/method.r)) / (Nτ + Nε))
        end
    end
    
    return efdr
end

"""
    calculate_efdr(method::PairedEFDR; show_progress=true)

Calculate paired empirical FDR with O(n²) complexity.
Uses progress monitoring that accounts for variable inner loop sizes.
"""
function calculate_efdr(method::PairedEFDR; show_progress::Bool=true)
    n = length(method.scores)
    efdr = zeros(eltype(method.qvals), n)
    
    # Validate inputs
    if !issorted(method.qvals)
        @warn "Q-values are not sorted. Results may be incorrect."
    end
    
    if any(method.complement_scores .< 0)
        @warn "Negative complement scores detected. Unpaired sequences?"
    end
    
    # Sort by q-value (ascending), then score (descending)
    sort_order = sortperm(collect(zip(method.qvals, -method.scores)))
    
    # Pre-compute total operations for accurate progress
    total_ops = sum(1:n)  # n*(n+1)/2
    completed_ops = 0
    
    # Pre-compute score relationships to avoid repeated comparisons
    println("Pre-computing score relationships...")
    score_rels = Vector{Symbol}(undef, n)
    
    for i in 1:n
        if method.is_original[i]
            score_rels[i] = :original
        else
            e_score = method.scores[i]
            o_score = method.complement_scores[i]
            
            if o_score < 0  # No complement
                score_rels[i] = :unpaired
            elseif e_score > o_score
                score_rels[i] = :entrap_wins
            elseif e_score == o_score
                score_rels[i] = :tie
            else
                score_rels[i] = :original_wins
            end
        end
    end
    
    # Main calculation with progress monitoring
    pb = show_progress ? ProgressBar(1:total_ops) : nothing
    if show_progress
        set_description(pb, "Calculating Paired EFDR")
    end
    
    for i in 1:n
        Nτ, Nε, Nεsτ, Nετs = 0, 0, 0, 0
        s = method.scores[sort_order[i]]
        
        for j in 1:i
            idx = sort_order[j]
            
            if method.is_original[idx]
                Nτ += 1
            else
                Nε += 1
                
                if score_rels[idx] == :entrap_wins
                    o_score = method.complement_scores[idx]
                    if o_score >= s
                        Nετs += 1
                    end
                elseif score_rels[idx] != :unpaired && method.scores[idx] >= s
                    o_score = method.complement_scores[idx]
                    if o_score < s
                        Nεsτ += 1
                    end
                end
            end
            
            # Update progress
            if show_progress
                completed_ops += 1
                update(pb)
            end
        end
        
        if Nτ + Nε > 0
            efdr[sort_order[i]] = min(1.0, (Nε + Nεsτ + 2*Nετs) / (Nτ + Nε))
        end
    end
    
    if show_progress
        finish!(pb)
    end
    
    return efdr
end
```

### 4.2 Monotonization
```julia
# src/core/scoring.jl

"""
    monotonize!(values::AbstractVector{T}) where T<:AbstractFloat

Ensure FDR values are monotonically non-decreasing.
Works in-place by traversing backwards through sorted results.

# Algorithm
Starting from the end (worst scores), ensures each FDR is at least
as large as the FDR of better-scoring results.
"""
function monotonize!(values::AbstractVector{T}) where T<:AbstractFloat
    current_min = T(1.0)
    
    for i in length(values):-1:1
        if values[i] > current_min
            values[i] = current_min
        else
            current_min = values[i]
        end
    end
    
    return values
end

# Handle missing values
function monotonize!(values::AbstractVector{Union{Missing, T}}) where T<:AbstractFloat
    current_min = T(1.0)
    
    for i in length(values):-1:1
        if ismissing(values[i])
            continue
        end
        
        if values[i] > current_min
            values[i] = current_min
        else
            current_min = values[i]
        end
    end
    
    return values
end
```

## Phase 5: Q-value Calculation

### 5.1 Q-value Calculations
```julia
# src/analysis/qvalue_calculation.jl

"""
    calculate_qvalues!(df::DataFrame; score_col=:PredVal)

Calculate local q-values using target-decoy approach.
Adds :local_qvalue column to the dataframe.
"""
function calculate_qvalues!(df::DataFrame; score_col::Symbol=:PredVal)
    println("Calculating local q-values...")
    
    # Validate columns
    if !hasproperty(df, score_col)
        error("Score column $score_col not found")
    end
    if !hasproperty(df, :decoy)
        error("Decoy column not found")
    end
    
    # Sort by score (descending), targets before decoys
    sort!(df, [score_col, :decoy], rev=[true, false])
    
    # Calculate q-values
    n = nrow(df)
    df[!, :local_qvalue] = zeros(Float32, n)
    
    Nτ, Nd = 0, 0
    
    for i in 1:n
        if df[i, :decoy]
            Nd += 1
        else
            Nτ += 1
        end
        
        if Nτ > 0
            df[i, :local_qvalue] = Nd / Nτ
        else
            df[i, :local_qvalue] = 0.0
        end
    end
    
    # Monotonize
    monotonize!(df[!, :local_qvalue])
    
    println("Local q-values: $(Nd) decoys, $(Nτ) targets")
    
    return nothing
end

"""
    calculate_global_qvalues!(df::DataFrame; score_col=:PredVal)

Calculate global q-values (best per precursor).
Adds :global_qvalue column to the dataframe.
"""
function calculate_global_qvalues!(df::DataFrame; score_col::Symbol=:PredVal)
    println("Calculating global q-values...")
    
    # Group by precursor
    gdf = groupby(df, [:decoy, :stripped_seq, :z])
    
    # Get best scoring per group
    global_summary = combine(gdf) do group
        idx = argmax(group[:, score_col])
        return group[idx:idx, [:decoy, :stripped_seq, :z, score_col]]
    end
    
    # Sort and calculate q-values
    sort!(global_summary, [score_col, :decoy], rev=[true, false])
    
    global_summary[!, :global_qvalue] = zeros(Float32, nrow(global_summary))
    
    Nτ, Nd = 0, 0
    
    for i in 1:nrow(global_summary)
        if global_summary[i, :decoy]
            Nd += 1
        else
            Nτ += 1
        end
        
        if Nτ > 0
            global_summary[i, :global_qvalue] = Nd / Nτ
        else
            global_summary[i, :global_qvalue] = 0.0
        end
    end
    
    monotonize!(global_summary[!, :global_qvalue])
    
    # Map back to original dataframe
    qval_dict = Dict(
        (row.decoy, row.stripped_seq, row.z) => row.global_qvalue
        for row in eachrow(global_summary)
    )
    
    df[!, :global_qvalue] = [qval_dict[(row.decoy, row.stripped_seq, row.z)] 
                              for row in eachrow(df)]
    
    println("Global q-values: $(Nd) decoys, $(Nτ) targets")
    
    return nothing
end
```

## Phase 6: Protein-Level Analysis (Detailed Implementation)

### 6.1 Protein Analysis Implementation
```julia
# src/analysis/protein_analysis.jl

"""
    prepare_protein_analysis(results_df::DataFrame, library_df::DataFrame; score_col=:PredVal)

Prepare data for protein-level analysis following the notebook approach:
1. Group by file_name, channel, decoy, entrapment_group, protein
2. Select best scoring PSM per protein group
3. Calculate protein-level q-values per run

This follows the exact logic from the notebook's protein analysis section.
"""
function prepare_protein_analysis(results_df::DataFrame, library_df::DataFrame; 
                                 score_col::Symbol=:PredVal)
    println("Preparing protein-level analysis...")
    
    # Step 1: Add entrapment_group column based on sequences
    target_seqs = Set(library_df[library_df.EntrapmentGroupId .== 0, :PeptideSequence])
    entrap_seqs = Set(library_df[library_df.EntrapmentGroupId .> 0, :PeptideSequence])
    
    results_df[!, :entrapment_group] = [seq ∈ entrap_seqs for seq in results_df.stripped_seq]
    
    # Step 2: Sort by PredVal (descending), targets above decoys for ties
    sort!(results_df, [score_col, :decoy], rev=[true, false])
    
    # Step 3: Group by run-channel-protein combination
    # Following notebook: [:file_name, :channel, :decoy, :entrapment_group, :protein]
    gresults_df = groupby(results_df, [:file_name, :channel, :decoy, :entrapment_group, :protein])
    
    # Step 4: Retain maximum scoring precursor for each protein group
    protein_df = combine(gresults_df) do group
        idx = argmax(group[:, score_col])
        return group[idx:idx, :]
    end
    
    # Step 5: Sort again for q-value calculation
    sort!(protein_df, [score_col, :decoy], rev=[true, false])
    
    # Step 6: Calculate protein q-values per run (exactly as notebook)
    protein_df[!, :Protein_Qvalue] = zeros(Float32, nrow(protein_df))
    
    for (run, results) in pairs(groupby(protein_df, :file_name))
        Nτ, Nd = 0, 0
        
        for i in 1:nrow(results)
            if results[i, :decoy]
                Nd += 1
            else
                Nτ += 1
            end
            results[i, :Protein_Qvalue] = Nd / Nτ
        end
        
        println("Run $(run.file_name): N_d $Nd N_t $Nτ")
        
        # Monotonize q-values
        monotonize!(results[!, :Protein_Qvalue])
    end
    
    return protein_df
end

"""
    calculate_protein_efdr!(protein_df::DataFrame, library_df::DataFrame, 
                           pair_dict::Dict, is_original_dict::Dict;
                           score_col=:PredVal)

Calculate entrapment FDR at the protein level using the same approach as precursors.
The key insight from the notebook is that protein-level EFDR uses:
- The rolled-up protein data (best PSM per protein)
- The same getPairedEmpiricalFdr function
- Grouping by file_name for per-run analysis
"""
function calculate_protein_efdr!(protein_df::DataFrame, library_df::DataFrame,
                                pair_dict::Dict, is_original_dict::Dict;
                                score_col::Symbol=:PredVal)
    
    # Remove decoys for EFDR calculation
    protein_no_decoys = filter(row -> !row.decoy, protein_df)
    
    # Initialize EFDR column
    protein_no_decoys[!, :protein_group_entrapment_fdr] = zeros(Float32, nrow(protein_no_decoys))
    
    # Group by file and calculate EFDR for each run
    gdf = groupby(protein_no_decoys, :file_name)
    
    for (key, run_df) in pairs(gdf)
        println("Calculating protein EFDR for run: $(key.file_name)")
        
        # Since we're using the same EFDR calculation as precursors,
        # we need to prepare the data in the same format
        
        # Get pairing information for proteins in this run
        n_proteins = nrow(run_df)
        is_original = Vector{Bool}(undef, n_proteins)
        pair_indices = Vector{Int}(undef, n_proteins)
        complement_scores = Vector{Float64}(undef, n_proteins)
        
        for i in 1:n_proteins
            seq = run_df[i, :stripped_seq]
            z = run_df[i, :z]
            key = (seq, z)
            
            is_original[i] = is_original_dict[key]
            pair_indices[i] = pair_dict[key]
            
            # Find complement score within this run
            complement_found = false
            for j in 1:n_proteins
                if i != j && pair_indices[j] == pair_indices[i] && is_original[i] != is_original[j]
                    complement_scores[i] = run_df[j, score_col]
                    complement_found = true
                    break
                end
            end
            
            if !complement_found
                complement_scores[i] = -1.0  # No complement in this run
            end
        end
        
        # Create PairedEFDR method for this run
        method = PairedEFDR(
            Float64.(run_df[!, score_col]),
            complement_scores,
            is_original,
            [is_original[i] ? 0 : 1 for i in 1:n_proteins],  # Simple entrap labels
            Float64.(run_df[!, :Protein_Qvalue]),
            pair_indices,
            1.0  # r_lib
        )
        
        # Calculate EFDR
        efdr_values = calculate_efdr(method; show_progress=false)
        monotonize!(efdr_values)
        
        # Store results
        run_df[!, :protein_group_entrapment_fdr] = efdr_values
    end
    
    return protein_no_decoys
end
```

## Phase 7: Main API Functions

### 7.1 Protein-Level Analysis API
```julia
# src/api.jl

"""
    run_protein_efdr_analysis(parquet_files, library_path; kwargs...)

Run empirical FDR analysis at the protein level following the notebook implementation.

# Process:
1. Load data and roll up to protein level
2. Calculate protein q-values per run
3. Calculate paired EFDR per run
4. Generate outputs with notebook-style plots

# Arguments
- `parquet_files`: Vector of paths to Parquet files with PSM results
- `library_path`: Path to TSV file with spectral library

# Keyword Arguments
- `output_dir`: Output directory (default: "efdr_output")
- `score_col`: Score column name (default: :PredVal)
- `r_lib`: Library to real entrapment ratio (default: 1.0)
- `show_progress`: Show progress bars (default: true)

# Returns
- DataFrame with protein-level EFDR analysis
"""
function run_protein_efdr_analysis(
    parquet_files::Vector{String},
    library_path::String;
    output_dir::String = "efdr_output",
    score_col::Symbol = :PredVal,
    r_lib::Float64 = 1.0,
    show_progress::Bool = true
)
    # Create output directory
    mkpath(output_dir)
    
    # Step 1: Load data
    results_df = load_parquet_results(parquet_files)
    library_df = load_spectral_library(library_path)
    
    # Step 2: Create pairing dictionaries (notebook style)
    println("Creating pairing dictionaries...")
    pair_dict = Dict{Tuple{String, Int}, Int}()
    is_original_dict = Dict{Tuple{String, Int}, Bool}()
    
    for i in 1:nrow(library_df)
        seq = library_df[i, :PeptideSequence]
        z = library_df[i, :PrecursorCharge]
        key = (seq, z)
        
        if !haskey(pair_dict, key)
            pair_dict[key] = library_df[i, :PrecursorIdx]
            is_original_dict[key] = library_df[i, :EntrapmentGroupId] == 0
        end
    end
    
    # Step 3: Prepare protein-level data
    protein_df = prepare_protein_analysis(results_df, library_df; score_col=score_col)
    
    # Step 4: Calculate protein EFDR
    protein_no_decoys = calculate_protein_efdr!(
        protein_df, library_df, pair_dict, is_original_dict;
        score_col=score_col
    )
    
    # Step 5: Generate outputs
    output_file = joinpath(output_dir, "protein_group_entrapment.tsv")
    CSV.write(output_file, protein_no_decoys, delim='\t')
    println("\nProtein results saved to: $output_file")
    
    # Step 6: Create visualization with notebook styling
    plot_protein_efdr_comparison(
        protein_no_decoys;
        output_path = joinpath(output_dir, "protein_efdr_comparison.pdf"),
        title = "Entrapment Analysis Protein Groups"
    )
    
    return protein_no_decoys
end

"""
    run_efdr_analysis(parquet_files, library_path; kwargs...)

Run empirical FDR analysis at the precursor level.
[Previous implementation remains the same, but ensure it uses notebook-style plotting]
"""
function run_efdr_analysis(
    parquet_files::Vector{String},
    library_path::String;
    output_dir::String = "efdr_output",
    score_col::Symbol = :PredVal,
    methods::Vector{DataType} = [CombinedEFDR, PairedEFDR],
    global_qval_threshold::Float64 = 0.01,
    r_lib::Float64 = 1.0,
    show_progress::Bool = true
)
    # [Previous implementation, but update plot call to use notebook style]
    # ...existing code...
    
    # Step 9: Create visualization with notebook styling
    plot_precursor_efdr_comparison(
        results_no_decoys;
        output_path = joinpath(output_dir, "precursor_efdr_comparison.pdf"),
        title = "Entrapment Analysis Precursors"
    )
    
    return results_no_decoys
end
```

## Phase 8: Visualization with Notebook Styling

### 8.1 Plotting Functions with Exact Notebook Style
```julia
# src/plotting/visualization.jl

"""
    plot_precursor_efdr_comparison(df; kwargs...)

Create EFDR comparison plot for precursors using exact notebook styling.
"""
function plot_precursor_efdr_comparison(
    df::DataFrame;
    output_path::String = "precursor_efdr_comparison.pdf",
    fdr_col::Symbol = :local_qvalue,
    efdr_col::Symbol = :precursor_entrapment_fdr,
    title::String = "Entrapment Analysis Precursors"
)
    # Exact plot styling from notebook
    p = plot(
        size = (400*1.5, 300*1.5),  # 600x450 pixels
        xlim = (0, 0.05),
        ylim = (0, 0.05),
        xlabel = "FDR",
        ylabel = "Entrapment FDR",
        title = title,
        titlefontsize = 16,
        xguidefontsize = 16,
        yguidefontsize = 16,
        tickfontsize = 12,
        legendfontsize = 12
    )
    
    # Diagonal reference line (exact style)
    plot!(p, [0, 0.05], [0, 0.05], 
          color = :black, 
          lw = 3, 
          label = nothing, 
          linestyle = :dash)
    
    # Plot each run with notebook's exact color
    for (key, run_df) in pairs(groupby(df, :file_name))
        plot!(p,
              run_df[!, fdr_col],
              run_df[!, efdr_col],
              lw = 3,
              label = nothing,  # No labels as per notebook
              color = RGB(0.39215686, 0.58431373, 0.92941176),  # Exact color
              alpha = 0.75)
    end
    
    savefig(p, output_path)
    println("Plot saved to: $output_path")
    
    return p
end

"""
    plot_protein_efdr_comparison(df; kwargs...)

Create EFDR comparison plot for protein groups using exact notebook styling.
"""
function plot_protein_efdr_comparison(
    df::DataFrame;
    output_path::String = "protein_efdr_comparison.pdf",
    fdr_col::Symbol = :Protein_Qvalue,
    efdr_col::Symbol = :protein_group_entrapment_fdr,
    title::String = "Entrapment Analysis Protein Groups"
)
    # Exact plot styling from notebook
    p = plot(
        size = (400*1.5, 300*1.5),  # 600x450 pixels
        xlim = (0, 0.05),
        ylim = (0, 0.05),
        xlabel = "FDR",
        ylabel = "Entrapment FDR",
        title = title
    )
    
    # Note: titlefontsize and other font sizes are set after the diagonal line in notebook
    plot!(p, [0, 0.05], [0, 0.05], 
          color = :black, 
          lw = 3, 
          label = nothing, 
          linestyle = :dash,
          titlefontsize = 16,
          xguidefontsize = 16,
          yguidefontsize = 16,
          tickfontsize = 12,
          legendfontsize = 12)
    
    # Plot each run
    for (key, run_df) in pairs(groupby(df, :file_name))
        plot!(p,
              run_df[!, fdr_col],
              run_df[!, efdr_col],
              lw = 3,
              label = nothing,
              color = RGB(0.39215686, 0.58431373, 0.92941176),
              alpha = 0.75)
    end
    
    savefig(p, output_path)
    println("Plot saved to: $output_path")
    
    return p
end
```

## Key Implementation Notes

### Protein-Level Analysis Details
The protein-level analysis follows these exact steps from the notebook:

1. **Grouping**: Group by `[:file_name, :channel, :decoy, :entrapment_group, :protein]`
2. **Best PSM Selection**: Take the highest scoring PSM per protein group
3. **Q-value Calculation**: Calculate per-run protein q-values (not global)
4. **EFDR Calculation**: Use the same paired EFDR algorithm, but on protein-rolled data
5. **Per-Run Analysis**: Always group by file_name for EFDR calculation

### Critical Differences from Module Approach
- Uses dictionary-based pairing for compatibility with notebook data
- Calculates protein q-values per run, not globally
- Maintains exact plot styling from notebook (size, colors, font sizes)
- No conversion to module's pairing system - uses notebook format throughout

### Plotting Style Notes
- Size: 600x450 pixels (400*1.5, 300*1.5)
- Color: RGB(0.39215686, 0.58431373, 0.92941176) - specific blue
- Alpha: 0.75 for transparency
- Line width: 3 for all lines
- No legend labels for individual runs
- Font sizes: title=16, axis labels=16, ticks=12

## Testing Considerations

### Test Data Requirements
- Must include protein column for protein-level testing
- Should have multiple runs (file_name values) to test per-run calculations
- Need paired sequences in both results and library

### Critical Test Cases
1. Protein rollup produces one row per protein/run/channel
2. Protein q-values are monotonic within each run
3. EFDR calculation handles missing pairs gracefully
4. Plot output matches notebook styling exactly

## Implementation Checklist (Updated)

- [ ] Remove all references to "TMT" from code and documentation
- [ ] Implement notebook-exact plot styling in visualization.jl
- [ ] Create protein_analysis.jl with detailed protein rollup logic
- [ ] Ensure protein EFDR uses same calculate_efdr method as precursors
- [ ] Test that protein q-values are calculated per-run, not globally
- [ ] Verify plot colors match RGB(0.39215686, 0.58431373, 0.92941176)
- [ ] Ensure font sizes match notebook (title=16, guides=16, ticks=12)
- [ ] Test with multi-run data to verify per-run calculations