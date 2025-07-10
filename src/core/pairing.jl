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
    compute_pairing_vectors!(library_df::DataFrame, results_df::DataFrame; kwargs...)

Add plex-aware pairing information directly to the results DataFrame.

Modifies `results_df` by adding the following columns:
- `:is_original` - Bool: true if peptide is original, false if entrapment
- `:pair_id` - Int: ID linking original/entrapment pairs from library
- `:entrap_label` - Int: 0 for original, 1 for entrapment (used in EFDR)
- `:complement_score` - Float32: plex & file-specific score of paired peptide (-1 if no pair)
- `:complement_indices` - Int: deprecated, always -1 (kept for compatibility)

The complement scores respect both file and plex boundaries, meaning the same
peptide in different plexes can have different complement scores.

# Arguments
- `library_df`: Library DataFrame with entrapment pairs
- `results_df`: Results DataFrame to modify (must have file_name and channel columns)

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
- Modified `results_df` with new columns added
"""
function compute_pairing_vectors!(
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
    
    # Step 3: Add remaining columns
    println("Adding entrap_label column (0=original, 1=entrapment)...")
    results_df[!, :entrap_label] = [orig ? 0 : 1 for orig in results_df[!, :is_original]]
    
    # Add complement_indices for backward compatibility
    results_df[!, :complement_indices] = fill(-1, nrow(results_df))
    
    println("Added pairing columns: is_original, pair_id, entrap_label, complement_score, complement_indices")
    
    return results_df
end

