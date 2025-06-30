# Implementation Plan: EntrapmentAnalyses.jl Module (Final Version)

## Reference Implementations

### Jupyter Notebook Reference
- **Path**: `/Users/nathanwamsley/Data/May2025/kmd_jmod_search/9plex/05202025/entrapment_test.ipynb`
- **Description**: Original implementation using dictionary-based pairing for 9-plex TMT data
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
- channel: Integer - TMT channel (optional)
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

## Phase 6: Protein-Level Analysis

### 6.1 Protein Analysis Implementation
```julia
# src/analysis/protein_analysis.jl

"""
    rollup_to_proteins!(df::DataFrame; score_col=:PredVal)

Roll up PSM results to protein groups.
Selects best scoring PSM per protein per run/channel.
"""
function rollup_to_proteins!(df::DataFrame; score_col::Symbol=:PredVal)
    println("Rolling up to protein groups...")
    
    # Validate columns
    required_cols = [:file_name, :channel, :decoy, :entrapment_group, :protein, score_col]
    missing_cols = setdiff(required_cols, names(df, Symbol))
    
    if !isempty(missing_cols)
        error("Missing required columns for protein rollup: $(join(missing_cols, ", "))")
    end
    
    # Group by run-channel-protein
    gdf = groupby(df, [:file_name, :channel, :decoy, :entrapment_group, :protein])
    
    # Get best PSM per protein group
    protein_df = combine(gdf) do group
        idx = argmax(group[:, score_col])
        return group[idx:idx, :]
    end
    
    println("Rolled up $(nrow(df)) PSMs to $(nrow(protein_df)) protein groups")
    
    return protein_df
end

"""
    calculate_protein_qvalues!(df::DataFrame; score_col=:PredVal)

Calculate protein-level q-values per run.
"""
function calculate_protein_qvalues!(df::DataFrame; score_col::Symbol=:PredVal)
    println("Calculating protein q-values...")
    
    # Sort for q-value calculation
    sort!(df, [score_col, :decoy], rev=[true, false])
    
    # Initialize column
    df[!, :protein_qvalue] = zeros(Float32, nrow(df))
    
    # Calculate per run
    for (key, run_df) in pairs(groupby(df, :file_name))
        Nτ, Nd = 0, 0
        
        for i in 1:nrow(run_df)
            if run_df[i, :decoy]
                Nd += 1
            else
                Nτ += 1
            end
            
            if Nτ > 0
                run_df[i, :protein_qvalue] = Nd / Nτ
            else
                run_df[i, :protein_qvalue] = 0.0
            end
        end
        
        monotonize!(run_df[!, :protein_qvalue])
        
        println("Run $(key.file_name): $(Nd) decoys, $(Nτ) targets")
    end
    
    return nothing
end
```

## Phase 7: Main API Functions

### 7.1 Precursor-Level Analysis
```julia
# src/api.jl

"""
    run_efdr_analysis(parquet_files, library_path; kwargs...)

Run empirical FDR analysis at the precursor level.

# Arguments
- `parquet_files`: Vector of paths to Parquet files with PSM results
- `library_path`: Path to TSV file with spectral library

# Keyword Arguments
- `output_dir`: Output directory (default: "efdr_output")
- `score_col`: Score column name (default: :PredVal)
- `methods`: EFDR methods to use (default: [CombinedEFDR, PairedEFDR])
- `global_qval_threshold`: Filter threshold (default: 0.01)
- `r_lib`: Library to real entrapment ratio (default: 1.0)
- `show_progress`: Show progress bars (default: true)

# Returns
- DataFrame with EFDR columns added
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
    # Create output directory
    mkpath(output_dir)
    
    # Step 1: Load data
    results_df = load_parquet_results(parquet_files)
    library_df = load_spectral_library(library_path)
    
    # Step 2: Add entrapment labels
    target_seqs = Set(library_df[library_df.EntrapmentGroupId .== 0, :PeptideSequence])
    entrap_seqs = Set(library_df[library_df.EntrapmentGroupId .> 0, :PeptideSequence])
    
    results_df[!, :entrapment_group] = [seq ∈ entrap_seqs for seq in results_df.stripped_seq]
    
    # Step 3: Calculate q-values
    calculate_qvalues!(results_df; score_col=score_col)
    calculate_global_qvalues!(results_df; score_col=score_col)
    
    # Step 4: Filter by global q-value
    n_before = nrow(results_df)
    filter!(row -> row.global_qvalue <= global_qval_threshold, results_df)
    n_after = nrow(results_df)
    println("Filtered to $(n_after)/$(n_before) PSMs at $(global_qval_threshold) global FDR")
    
    # Step 5: Remove decoys for EFDR calculation
    results_no_decoys = filter(row -> !row.decoy, results_df)
    println("Analyzing $(nrow(results_no_decoys)) target PSMs")
    
    # Step 6: Compute pairing vectors
    pairing_info = compute_pairing_vectors(library_df, results_no_decoys; 
                                          show_progress=show_progress)
    
    # Step 7: Group by run and calculate EFDR
    results_no_decoys[!, :precursor_entrapment_fdr] = zeros(Float32, nrow(results_no_decoys))
    
    for (key, run_df) in pairs(groupby(results_no_decoys, :file_name))
        println("\nProcessing run: $(key.file_name)")
        
        # Get indices for this run
        run_indices = run_df.index
        
        # Extract data for this run
        scores = Float64.(run_df[!, score_col])
        qvals = Float64.(run_df[!, :local_qvalue])
        
        # Get complement scores
        complement_scores = Vector{Float64}(undef, length(run_indices))
        for (i, idx) in enumerate(run_indices)
            comp_idx = pairing_info.complement_indices[idx]
            if comp_idx > 0
                complement_scores[i] = results_no_decoys[comp_idx, score_col]
            else
                complement_scores[i] = -1.0  # Indicates no pair
            end
        end
        
        # Calculate paired EFDR
        method = PairedEFDR(
            scores,
            complement_scores,
            pairing_info.is_original[run_indices],
            pairing_info.entrap_labels[run_indices],
            qvals,
            pairing_info.pair_indices[run_indices],
            r_lib
        )
        
        efdr_values = calculate_efdr(method; show_progress=show_progress)
        monotonize!(efdr_values)
        
        # Store results
        run_df[!, :precursor_entrapment_fdr] = efdr_values
    end
    
    # Step 8: Generate outputs
    output_file = joinpath(output_dir, "precursor_entrapment_results.tsv")
    CSV.write(output_file, results_no_decoys, delim='\t')
    println("\nResults saved to: $output_file")
    
    # Step 9: Create visualization
    plot_efdr_comparison(
        results_no_decoys;
        output_path = joinpath(output_dir, "precursor_efdr_comparison.pdf")
    )
    
    return results_no_decoys
end

"""
    run_protein_efdr_analysis(parquet_files, library_path; kwargs...)

Run empirical FDR analysis at the protein level.
"""
function run_protein_efdr_analysis(
    parquet_files::Vector{String},
    library_path::String;
    output_dir::String = "efdr_output",
    score_col::Symbol = :PredVal,
    r_lib::Float64 = 1.0,
    show_progress::Bool = true
)
    # Load and process as before
    results_df = load_parquet_results(parquet_files)
    library_df = load_spectral_library(library_path)
    
    # Add entrapment labels
    entrap_seqs = Set(library_df[library_df.EntrapmentGroupId .> 0, :PeptideSequence])
    results_df[!, :entrapment_group] = [seq ∈ entrap_seqs for seq in results_df.stripped_seq]
    
    # Roll up to proteins
    protein_df = rollup_to_proteins!(results_df; score_col=score_col)
    
    # Calculate protein q-values
    calculate_protein_qvalues!(protein_df; score_col=score_col)
    
    # Remove decoys
    protein_no_decoys = filter(row -> !row.decoy, protein_df)
    
    # Compute pairing and EFDR (similar process as precursor level)
    # ... (implementation follows same pattern)
    
    return protein_no_decoys
end
```

## Phase 8: Visualization

### 8.1 Plotting Functions
```julia
# src/plotting/visualization.jl

"""
    plot_efdr_comparison(df; kwargs...)

Create EFDR comparison plot.
"""
function plot_efdr_comparison(
    df::DataFrame;
    output_path::String = "efdr_comparison.pdf",
    fdr_col::Symbol = :local_qvalue,
    efdr_col::Symbol = :precursor_entrapment_fdr,
    title::String = "Entrapment FDR Analysis"
)
    p = plot(
        size = (600, 450),
        xlim = (0, 0.05),
        ylim = (0, 0.05),
        xlabel = "FDR",
        ylabel = "Entrapment FDR",
        title = title,
        titlefontsize = 16,
        guidefontsize = 16,
        tickfontsize = 12,
        legendfontsize = 12
    )
    
    # Diagonal reference line
    plot!(p, [0, 0.05], [0, 0.05], 
          color = :black, 
          lw = 3, 
          label = nothing, 
          linestyle = :dash)
    
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

## Phase 9: Testing Infrastructure

### 9.1 Unit Tests
```julia
# test/unit/test_pairing.jl

using Test
using EntrapmentAnalyses

@testset "Pairing Tests" begin
    # Create test data
    library_df = DataFrame(
        PeptideSequence = ["PEPTIDE", "EPPTIDE", "SEQUENCE", "SEQEUNCE"],
        PrecursorCharge = [2, 2, 3, 3],
        EntrapmentGroupId = [0, 1, 0, 1],
        PrecursorIdx = [1, 1, 2, 2]
    )
    
    results_df = DataFrame(
        stripped_seq = ["PEPTIDE", "EPPTIDE", "SEQUENCE"],
        z = [2, 2, 3]
    )
    
    @testset "compute_pairing_vectors" begin
        pairing = compute_pairing_vectors(library_df, results_df; 
                                        show_progress=false)
        
        @test pairing.is_original == [true, false, true]
        @test pairing.pair_indices == [1, 1, 2]
        @test pairing.complement_indices == [2, 1, -1]  # Last has no pair in results
    end
end

# test/unit/test_efdr_methods.jl

@testset "EFDR Methods" begin
    @testset "monotonize!" begin
        values = [0.01, 0.005, 0.02, 0.015, 0.03]
        monotonize!(values)
        @test issorted(values)
        @test values == [0.005, 0.005, 0.015, 0.015, 0.03]
    end
    
    @testset "CombinedEFDR" begin
        scores = [0.9, 0.8, 0.7, 0.6, 0.5]
        entrap_labels = [0, 1, 0, 1, 0]
        qvals = [0.001, 0.002, 0.003, 0.004, 0.005]
        
        method = CombinedEFDR(scores, entrap_labels, qvals, 1.0)
        efdr = calculate_efdr(method; show_progress=false)
        
        @test length(efdr) == length(scores)
        @test all(0 .<= efdr .<= 1)
    end
end
```

### 9.2 Integration Tests
```julia
# test/integration/test_full_pipeline.jl

@testset "Full Pipeline Integration" begin
    # Use small test files
    test_dir = joinpath(@__DIR__, "test_data")
    parquet_files = [joinpath(test_dir, "sample_results.parquet")]
    library_path = joinpath(test_dir, "sample_library.tsv")
    
    @testset "Precursor Analysis" begin
        results = run_efdr_analysis(
            parquet_files,
            library_path;
            output_dir = mktempdir(),
            show_progress = false
        )
        
        @test hasproperty(results, :precursor_entrapment_fdr)
        @test all(0 .<= results.precursor_entrapment_fdr .<= 1)
        @test issorted(results.precursor_entrapment_fdr)
    end
end
```

### 9.3 Test Runner
```julia
# test/runtests.jl

using Test
using EntrapmentAnalyses

# Unit tests
include("unit/test_pairing.jl")
include("unit/test_efdr_methods.jl")
include("unit/test_qvalue.jl")

# Integration tests
include("integration/test_full_pipeline.jl")
```

## Implementation Checklist

- [ ] Phase 1: Core module setup
  - [ ] Create Project.toml with dependencies
  - [ ] Set up module structure
  - [ ] Create directory structure

- [ ] Phase 2: Data loading
  - [ ] Implement Parquet loader with validation
  - [ ] Implement TSV library loader
  - [ ] Add error handling for missing files/columns

- [ ] Phase 3: Pairing system
  - [ ] Implement compute_pairing_vectors
  - [ ] Add progress monitoring
  - [ ] Validate pairing completeness

- [ ] Phase 4: EFDR methods
  - [ ] Implement CombinedEFDR struct and calculation
  - [ ] Implement PairedEFDR with progress monitoring
  - [ ] Add monotonize! function

- [ ] Phase 5: Q-value calculations
  - [ ] Implement local q-value calculation
  - [ ] Implement global q-value calculation
  - [ ] Add filtering capabilities

- [ ] Phase 6: Protein analysis
  - [ ] Implement protein rollup
  - [ ] Add protein q-value calculation
  - [ ] Integrate with EFDR pipeline

- [ ] Phase 7: API functions
  - [ ] Implement run_efdr_analysis
  - [ ] Implement run_protein_efdr_analysis
  - [ ] Add output generation

- [ ] Phase 8: Visualization
  - [ ] Implement EFDR comparison plots
  - [ ] Add multi-run support
  - [ ] Match notebook styling

- [ ] Phase 9: Testing
  - [ ] Create unit tests
  - [ ] Create integration tests
  - [ ] Generate test data

## Success Metrics

1. **Correctness**: Results match notebook implementation within numerical precision
2. **Performance**: Paired EFDR calculation completes in reasonable time with progress
3. **Usability**: Clear error messages and progress indicators
4. **Compatibility**: Accepts exact same input format as notebook
5. **Testing**: >90% code coverage with comprehensive tests