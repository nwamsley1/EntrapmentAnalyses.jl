# Implementation Plan: Plex-Specific Original Score Calculation

## Overview

This plan updates the pairing system to calculate plex-specific complement scores while keeping the EFDR calculation methods unchanged. The key insight is that pairing dictionaries are built once from the library, but the plex score dictionary is built per file.

## Key Architecture

1. **Global** (built once): `pair_dict` and `is_original_dict` from library
2. **Per-file**: `plex_prec_to_scores` dictionary
3. **Result**: Add `complement_score` column to the dataframe

## Core Changes

### 1. New Type Definitions Using Structs

```julia
# Type for peptide identification
struct PeptideKey
    mod_seq::String
    z::UInt8
end

# Type for plex-specific pairing within a file
struct PlexPairKey
    plex::UInt8
    pair_id::Int64
end

# Type for score pairs
struct ScorePair
    original_score::Float32
    entrapment_score::Float32
end

# Make types hashable and comparable
Base.hash(k::PeptideKey, h::UInt) = hash((k.mod_seq, k.z), h)
Base.:(==)(a::PeptideKey, b::PeptideKey) = a.mod_seq == b.mod_seq && a.z == b.z

Base.hash(k::PlexPairKey, h::UInt) = hash((k.plex, k.pair_id), h)
Base.:(==)(a::PlexPairKey, b::PlexPairKey) = a.plex == b.plex && a.pair_id == b.pair_id
```

### 2. New Functions in `pairing.jl`

#### NEW: `init_entrapment_pairs_dict` (called once globally)

```julia
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
```

#### NEW: `add_plex_complement_scores!` (modifies dataframe in place)

```julia
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
        results_df[i, :is_original] = is_original_dict[peptide_key]
        results_df[i, :pair_id] = pair_dict[peptide_key]
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
```

#### REPLACE: `compute_pairing_vectors` with new version

```julia
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
```

### 3. Alternative Approach: Process During EFDR Calculation

If we want to avoid modifying the global dataframe, we can build the plex scores on-the-fly:

```julia
function get_plex_complement_scores_for_run(
    run_df::DataFrame,
    pair_dict::Dict{PeptideKey, Int64},
    is_original_dict::Dict{PeptideKey, Bool};
    score_col::Symbol = :PredVal,
    channel_col::Symbol = :channel,
    seq_col::Symbol = :stripped_seq,
    charge_col::Symbol = :z
)
    # Build plex score dictionary for this specific run
    plex_prec_to_scores = Dict{PlexPairKey, ScorePair}()
    
    # First pass: build scores
    for i in 1:nrow(run_df)
        peptide_key = PeptideKey(run_df[i, seq_col], UInt8(run_df[i, charge_col]))
        plex = UInt8(run_df[i, channel_col])
        pair_id = pair_dict[peptide_key]
        plex_key = PlexPairKey(plex, pair_id)
        is_original = is_original_dict[peptide_key]
        
        if !haskey(plex_prec_to_scores, plex_key)
            plex_prec_to_scores[plex_key] = ScorePair(-1.0f0, -1.0f0)
        end
        
        scores = plex_prec_to_scores[plex_key]
        
        if is_original
            plex_prec_to_scores[plex_key] = ScorePair(Float32(run_df[i, score_col]), scores.entrapment_score)
        else
            plex_prec_to_scores[plex_key] = ScorePair(scores.original_score, Float32(run_df[i, score_col]))
        end
    end
    
    # Second pass: extract complement scores
    complement_scores = Vector{Float32}(undef, nrow(run_df))
    
    for i in 1:nrow(run_df)
        peptide_key = PeptideKey(run_df[i, seq_col], UInt8(run_df[i, charge_col]))
        plex = UInt8(run_df[i, channel_col])
        pair_id = pair_dict[peptide_key]
        plex_key = PlexPairKey(plex, pair_id)
        is_original = is_original_dict[peptide_key]
        
        if haskey(plex_prec_to_scores, plex_key)
            scores = plex_prec_to_scores[plex_key]
            complement_scores[i] = is_original ? scores.entrapment_score : scores.original_score
        else
            complement_scores[i] = -1.0f0
        end
    end
    
    return complement_scores
end
```

### 4. Usage in API

#### Option A: Modify dataframe once (recommended)

```julia
# In run_efdr_analysis:
# Build pairing info once for entire dataset
pairing_info = compute_pairing_vectors(library_df, results_no_decoys; show_progress=show_progress)

# Now results_no_decoys has complement_score column with file & plex specific values
# When processing by file, the complement scores are already correct
for (key, run_df) in pairs(groupby(results_no_decoys, :file_name))
    # complement_scores are already in the dataframe
    paired_efdr_values = calculate_paired_efdr(
        Float64.(run_df[!, score_col]),
        Float64.(run_df[!, :complement_score]),  # Already file & plex specific!
        run_df[!, :is_original],
        Float64.(run_df[!, :local_qvalue]);
        r = r_lib,
        show_progress = show_progress
    )
end
```

#### Option B: Calculate on-the-fly

```julia
# Build global dictionaries once
pair_dict, is_original_dict = init_entrapment_pairs_dict(library_df)

# Process each file
for (key, run_df) in pairs(groupby(results_no_decoys, :file_name))
    # Get complement scores for this run
    complement_scores = get_plex_complement_scores_for_run(
        run_df, pair_dict, is_original_dict
    )
    
    # Calculate EFDR
    paired_efdr_values = calculate_paired_efdr(
        Float64.(run_df[!, score_col]),
        Float64.(complement_scores),
        # ... other params
    )
end
```

## Benefits of This Approach

1. **Efficiency**: Library dictionaries built once, not per file
2. **Correctness**: Plex scores are file-specific as in the notebook
3. **Flexibility**: Can either modify dataframe or calculate on-the-fly
4. **Clear separation**: Global pairing info vs file-specific scores

## What Stays the Same

- All EFDR calculation methods remain unchanged
- The interface mostly stays the same
- Backward compatibility is maintained

## Summary

The key insight is separating:
1. **Global information** from library (pair IDs, original/entrapment status)
2. **File-specific information** (actual scores within each file/plex combination)

This matches the notebook's behavior while being more efficient than rebuilding dictionaries for each file.