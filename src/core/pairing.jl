# Efficient pairing system for entrapment analysis

using DataFrames
using ProgressBars

# Type definitions for plex-specific pairing

"""
Type for peptide identification based on sequence and charge.
"""
struct PeptideKey
    mod_seq::String
    z::UInt8
end

"""
Type for plex-specific pairing within a file.
"""
struct PlexPairKey
    plex::UInt8
    pair_id::Int64
end

"""
Type for storing original and entrapment scores.
"""
struct ScorePair
    original_score::Float32
    entrapment_score::Float32
end

# Make types hashable and comparable
Base.hash(k::PeptideKey, h::UInt) = hash((k.mod_seq, k.z), h)
Base.:(==)(a::PeptideKey, b::PeptideKey) = a.mod_seq == b.mod_seq && a.z == b.z

Base.hash(k::PlexPairKey, h::UInt) = hash((k.plex, k.pair_id), h)
Base.:(==)(a::PlexPairKey, b::PlexPairKey) = a.plex == b.plex && a.pair_id == b.pair_id

"""
    compute_pairing_vectors_old(lib_sequences, lib_charges, lib_entrap_groups, 
                           lib_pair_ids, results_sequences, results_charges; 
                           show_progress=true)

[DEPRECATED] Old version of compute_pairing_vectors without plex-specific pairing.
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
function compute_pairing_vectors_old(
    lib_sequences::AbstractVector{<:AbstractString},
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
    compute_pairing_vectors_old(library_df::DataFrame, results_df::DataFrame; kwargs...)

[DEPRECATED] DataFrame wrapper for old compute_pairing_vectors that extracts the appropriate columns.

# Keyword Arguments
- `lib_seq_col`: Library sequence column (default: :PeptideSequence)
- `lib_charge_col`: Library charge column (default: :PrecursorCharge)
- `lib_entrap_col`: Library entrapment group column (default: :EntrapmentGroupId)
- `lib_pair_col`: Library pair ID column (default: :PrecursorIdx)
- `results_seq_col`: Results sequence column (default: :stripped_seq)
- `results_charge_col`: Results charge column (default: :z)
- `show_progress`: Show progress bars (default: true)
"""
function compute_pairing_vectors_old(
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
    return compute_pairing_vectors_old(
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

[DEPRECATED] Extract complement scores based on complement indices.
Returns -1.0 for entries with no complement (complement_indices[i] == -1).

This function is deprecated. Use compute_pairing_vectors which now returns
plex-specific complement_scores in the result tuple.
"""
function get_complement_scores(scores::AbstractVector{T}, complement_indices::AbstractVector{Int}) where T<:Real
    @warn "get_complement_scores is deprecated. Use compute_pairing_vectors which returns plex-specific complement_scores." maxlog=1
    
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

# New functions for plex-specific pairing

"""
    init_entrapment_pairs_dict(lib_df::DataFrame; kwargs...)

Initialize dictionaries from the spectral library for pairing information.

# Arguments
- `lib_df`: Library DataFrame

# Keyword Arguments
- `mod_seq_col`: Column for peptide sequence (default: :PeptideSequence)
- `charge_col`: Column for charge state (default: :PrecursorCharge)
- `entrap_group_col`: Column for entrapment group ID (default: :EntrapmentGroupId)
- `pair_column`: Column for pair ID (default: :PrecursorIdx)

# Returns
Tuple of (pair_dict, is_original_dict)
"""
function init_entrapment_pairs_dict(
    lib_df::DataFrame;
    mod_seq_col::Symbol = :PeptideSequence,
    charge_col::Symbol = :PrecursorCharge,
    entrap_group_col::Symbol = :EntrapmentGroupId,
    pair_column::Symbol = :PrecursorIdx
)
    # Create dictionary for pair IDs
    pair_dict = Dict{PeptideKey, Int64}()
    
    # Create dictionary for original/entrapment status
    is_original_dict = Dict{PeptideKey, Bool}()
    
    for row_idx in 1:nrow(lib_df)
        mod_seq_value = lib_df[row_idx, mod_seq_col]  
        z = UInt8(lib_df[row_idx, charge_col])
        key = PeptideKey(mod_seq_value, z)
        
        if !haskey(pair_dict, key)
            pair_dict[key] = lib_df[row_idx, pair_column]
            is_original_dict[key] = iszero(lib_df[row_idx, entrap_group_col])
        end
    end
    
    return pair_dict, is_original_dict
end

"""
    add_plex_complement_scores!(results_df::DataFrame, pair_dict, is_original_dict; kwargs...)

Add plex-specific complement scores to the results DataFrame.

Modifies the DataFrame in place by adding:
- `complement_score`: Plex-specific complement scores (-1 if no pair)
- `is_original`: Boolean indicating if peptide is original
- `pair_id`: Pair identifier from library

# Arguments
- `results_df`: Results DataFrame to modify
- `pair_dict`: Dictionary mapping PeptideKey to pair IDs
- `is_original_dict`: Dictionary mapping PeptideKey to original status

# Keyword Arguments
- `score_col`: Score column (default: :PredVal)
- `channel_col`: Channel/plex column (default: :channel)
- `seq_col`: Sequence column (default: :stripped_seq)
- `charge_col`: Charge column (default: :z)
- `file_col`: File name column (default: :file_name)
- `show_progress`: Show progress (default: true)
"""
function add_plex_complement_scores!(
    results_df::DataFrame,
    pair_dict::Dict{PeptideKey, Int64},
    is_original_dict::Dict{PeptideKey, Bool};
    score_col::Symbol = :PredVal,
    channel_col::Symbol = :channel,
    seq_col::Symbol = :stripped_seq,
    charge_col::Symbol = :z,
    file_col::Symbol = :file_name,
    show_progress::Bool = true
)
    # Add complement score column if it doesn't exist
    if !hasproperty(results_df, :complement_score)
        results_df[!, :complement_score] = fill(-1.0f0, nrow(results_df))
    end
    
    # Add is_original and pair_id columns for later use
    if !hasproperty(results_df, :is_original)
        results_df[!, :is_original] = Vector{Bool}(undef, nrow(results_df))
    end
    if !hasproperty(results_df, :pair_id)
        results_df[!, :pair_id] = Vector{Int}(undef, nrow(results_df))
    end
    
    # First pass: populate is_original and pair_id for all rows
    for i in 1:nrow(results_df)
        peptide_key = PeptideKey(results_df[i, seq_col], UInt8(results_df[i, charge_col]))
        if haskey(pair_dict, peptide_key)
            results_df[i, :is_original] = is_original_dict[peptide_key]
            results_df[i, :pair_id] = pair_dict[peptide_key]
        else
            error("Sequence not found in library: $(results_df[i, seq_col]) with charge $(results_df[i, charge_col])")
        end
    end
    
    # Process each file separately
    file_groups = groupby(results_df, file_col)
    
    for (file_key, file_df) in pairs(file_groups)
        if show_progress
            println("Processing file: $(file_key[file_col])")
        end
        
        # Build plex score dictionary for this file
        plex_prec_to_scores = Dict{PlexPairKey, ScorePair}()
        
        # First pass: build the score dictionary
        for row in eachrow(file_df)
            plex = UInt8(row[channel_col])
            pair_id = row[:pair_id]
            plex_key = PlexPairKey(plex, pair_id)
            
            # Initialize if needed
            if !haskey(plex_prec_to_scores, plex_key)
                plex_prec_to_scores[plex_key] = ScorePair(-1.0f0, -1.0f0)
            end
            
            # Get current scores
            scores = plex_prec_to_scores[plex_key]
            
            # Update the appropriate score
            if row[:is_original]
                plex_prec_to_scores[plex_key] = ScorePair(
                    Float32(row[score_col]), 
                    scores.entrapment_score
                )
            else
                plex_prec_to_scores[plex_key] = ScorePair(
                    scores.original_score,
                    Float32(row[score_col])
                )
            end
        end
        
        # Second pass: assign complement scores
        for row in eachrow(file_df)
            plex = UInt8(row[channel_col])
            pair_id = row[:pair_id]
            plex_key = PlexPairKey(plex, pair_id)
            
            if haskey(plex_prec_to_scores, plex_key)
                scores = plex_prec_to_scores[plex_key]
                if row[:is_original]
                    row[:complement_score] = scores.entrapment_score
                else
                    row[:complement_score] = scores.original_score
                end
            end
            # If not found, remains -1
        end
    end
    
    return results_df
end

"""
    compute_pairing_vectors(library_df::DataFrame, results_df::DataFrame; kwargs...)

Compute plex-aware pairing vectors with file and plex-specific complement scores.

# Arguments
- `library_df`: Library DataFrame with entrapment pairs
- `results_df`: Results DataFrame with PSM data

# Keyword Arguments
- `lib_seq_col`: Library sequence column (default: :PeptideSequence)
- `lib_charge_col`: Library charge column (default: :PrecursorCharge)
- `lib_entrap_col`: Library entrapment group column (default: :EntrapmentGroupId)
- `lib_pair_col`: Library pair ID column (default: :PrecursorIdx)
- `results_seq_col`: Results sequence column (default: :stripped_seq)
- `results_charge_col`: Results charge column (default: :z)
- `channel_col`: Channel/plex column (default: :channel)
- `score_col`: Score column (default: :PredVal)
- `file_col`: File name column (default: :file_name)
- `show_progress`: Show progress bars (default: true)

# Returns
Named tuple with:
- `is_original`: Bool vector indicating if each result is an original sequence
- `pair_indices`: Int vector with pair IDs for each result
- `entrap_labels`: Int vector with entrapment group IDs
- `complement_indices`: Int vector (deprecated, always -1)
- `complement_scores`: Float32 vector with plex-specific complement scores
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
    channel_col::Symbol = :channel,
    score_col::Symbol = :PredVal,
    file_col::Symbol = :file_name,
    show_progress::Bool = true
)
    println("Computing plex-aware pairing vectors...")
    
    # Step 1: Initialize dictionaries from library (done once)
    pair_dict, is_original_dict = init_entrapment_pairs_dict(
        library_df;
        mod_seq_col = lib_seq_col,
        charge_col = lib_charge_col,
        entrap_group_col = lib_entrap_col,
        pair_column = lib_pair_col
    )
    
    # Step 2: Add complement scores to the dataframe (handles per-file logic internally)
    add_plex_complement_scores!(
        results_df,
        pair_dict,
        is_original_dict;
        score_col = score_col,
        channel_col = channel_col,
        seq_col = results_seq_col,
        charge_col = results_charge_col,
        file_col = file_col,
        show_progress = show_progress
    )
    
    # Step 3: Extract vectors for compatibility with existing code
    n_results = nrow(results_df)
    is_original_vec = results_df[!, :is_original]
    pair_indices_vec = results_df[!, :pair_id]
    complement_scores_vec = results_df[!, :complement_score]
    
    # Create entrap_labels from is_original
    entrap_labels_vec = [orig ? 0 : 1 for orig in is_original_vec]
    
    # For backward compatibility
    complement_indices = fill(-1, n_results)
    
    return (
        is_original = is_original_vec,
        pair_indices = pair_indices_vec,
        entrap_labels = entrap_labels_vec,
        complement_indices = complement_indices,  # Deprecated
        complement_scores = complement_scores_vec # File & plex specific!
    )
end