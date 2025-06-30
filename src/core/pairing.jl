# Efficient pairing system for entrapment analysis

"""
    compute_pairing_vectors(lib_sequences, lib_charges, lib_entrap_groups, 
                           lib_pair_ids, results_sequences, results_charges; 
                           show_progress=true)

Core function to compute pairing information as vectors for efficient EFDR calculation.

# Arguments
- `lib_sequences`: Library peptide sequences
- `lib_charges`: Library charge states
- `lib_entrap_groups`: Entrapment group IDs (0 = original, >0 = entrapment)
- `lib_pair_ids`: Pair identifiers linking originals to entrapments
- `results_sequences`: Result peptide sequences
- `results_charges`: Result charge states

# Returns
Named tuple with:
- `is_original`: Bool vector indicating if each result is an original sequence
- `pair_indices`: Int vector with pair IDs for each result
- `entrap_labels`: Int vector with entrapment group IDs
- `complement_indices`: Int vector with indices of paired sequences (-1 if no pair)
"""
function compute_pairing_vectors(
    lib_sequences::AbstractVector{String},
    lib_charges::AbstractVector{<:Integer},
    lib_entrap_groups::AbstractVector{<:Integer},
    lib_pair_ids::AbstractVector{<:Integer},
    results_sequences::AbstractVector{String},
    results_charges::AbstractVector{<:Real};
    show_progress::Bool = true
)
    n_results = length(results_sequences)
    n_lib = length(lib_sequences)
    
    # Validate input lengths
    if length(lib_charges) != n_lib || length(lib_entrap_groups) != n_lib || length(lib_pair_ids) != n_lib
        error("Library vectors must have the same length")
    end
    if length(results_charges) != n_results
        error("Result sequence and charge vectors must have the same length")
    end
    
    println("Computing pairing vectors for $n_results PSMs...")
    
    # Pre-allocate output vectors
    is_original = Vector{Bool}(undef, n_results)
    pair_indices = Vector{Int}(undef, n_results)
    entrap_labels = Vector{Int}(undef, n_results)
    complement_indices = fill(-1, n_results)
    
    # Build efficient lookup structure
    lib_lookup = Dict{Tuple{String, Int}, Int}()
    for i in 1:n_lib
        key = (lib_sequences[i], lib_charges[i])
        lib_lookup[key] = i
    end
    
    # Fill basic pairing information
    pb = show_progress ? ProgressBar(1:n_results) : (1:n_results)
    if show_progress
        set_description(pb, "Mapping sequences to library")
    end
    
    for i in pb
        key = (results_sequences[i], results_charges[i])
        
        if haskey(lib_lookup, key)
            lib_idx = lib_lookup[key]
            is_original[i] = lib_entrap_groups[lib_idx] == 0
            pair_indices[i] = lib_pair_ids[lib_idx]
            entrap_labels[i] = lib_entrap_groups[lib_idx]
        else
            error("Sequence not found in library: $(results_sequences[i]) with charge $(results_charges[i])")
        end
    end
    
    # Build pair lookup for complement finding
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
    pb2 = show_progress ? ProgressBar(1:n_results) : (1:n_results)
    if show_progress
        set_description(pb2, "Pairing complements")
    end
    
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

"""
    compute_pairing_vectors(library_df::DataFrame, results_df::DataFrame; kwargs...)

DataFrame wrapper for compute_pairing_vectors that extracts the appropriate columns.

# Keyword Arguments
- `lib_seq_col`: Library sequence column (default: :PeptideSequence)
- `lib_charge_col`: Library charge column (default: :PrecursorCharge)
- `lib_entrap_col`: Library entrapment group column (default: :EntrapmentGroupId)
- `lib_pair_col`: Library pair ID column (default: :PrecursorIdx)
- `results_seq_col`: Results sequence column (default: :stripped_seq)
- `results_charge_col`: Results charge column (default: :z)
- `show_progress`: Show progress bars (default: true)
"""
function compute_pairing_vectors(
    library_df::DataFrame,
    results_df::DataFrame;
    lib_seq_col::Symbol = :PeptideSequence,
    lib_charge_col::Symbol = :PrecursorCharge,
    lib_entrap_col::Symbol = :EntrapmentGroupId,
    lib_pair_col::Symbol = :PrecursorIdx,
    results_seq_col::Symbol = :stripped_seq,
    results_charge_col::Symbol = :z,
    show_progress::Bool = true
)
    return compute_pairing_vectors(
        library_df[!, lib_seq_col],
        library_df[!, lib_charge_col],
        library_df[!, lib_entrap_col],
        library_df[!, lib_pair_col],
        results_df[!, results_seq_col],
        results_df[!, results_charge_col];
        show_progress = show_progress
    )
end

"""
    get_complement_scores(scores::AbstractVector{T}, complement_indices::AbstractVector{Int}) where T

Extract complement scores based on complement indices.
Returns -1.0 for entries with no complement (complement_indices[i] == -1).
"""
function get_complement_scores(scores::AbstractVector{T}, complement_indices::AbstractVector{Int}) where T<:Real
    n = length(scores)
    if length(complement_indices) != n
        error("Scores and complement_indices must have the same length")
    end
    
    complement_scores = Vector{T}(undef, n)
    
    for i in 1:n
        if complement_indices[i] > 0 && complement_indices[i] <= n
            complement_scores[i] = scores[complement_indices[i]]
        else
            complement_scores[i] = T(-1)  # No complement
        end
    end
    
    return complement_scores
end