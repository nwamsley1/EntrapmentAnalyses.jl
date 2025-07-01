# PlexId/Channel Pairing Analysis: Notebook vs Repository

## Executive Summary

The notebook implementation uses **plex-specific pairing** for EFDR calculations, while the repository implementation uses **global pairing** without considering plex/channel information. This is a critical difference that could lead to incorrect EFDR calculations for multiplex experiments.

## Detailed Analysis

### 1. Notebook Implementation (plex_prec_to_scores)

In the notebook's `getPairedEmpiricalFdr` function, the pairing is plex-specific:

```julia
# Key structure includes BOTH plex and pair_id
plex_prec_to_scores = Dictionary{
    @NamedTuple{plex::UInt8, pair_id::Int64},
    @NamedTuple{original_score::Float32, entrapment_score::Float32}
}()

# Later in the code:
plex = UInt8(df[i,channel_id_col])
pair_id = pair_id_dict[(mod_seq = mod_seq, z = z)]
key = (plex = plex, pair_id = pair_id)  # Plex-specific key!
```

**Key Points:**
- Each plex (multiplex channel) maintains separate score pairs
- The same peptide in different plexes are treated as separate entities
- This allows for plex-specific competition between original and entrapment peptides

### 2. Repository Implementation (compute_pairing_vectors)

In the repository's implementation, pairing is global across all plexes:

```julia
# From compute_pairing_vectors in pairing.jl
lib_lookup = Dict{Tuple{String, Int}, Int}()
for i in 1:n_lib
    key = (lib_sequences[i], lib_charges[i])  # No plex/channel!
    lib_lookup[key] = i
end
```

**Key Points:**
- Pairing is based only on (sequence, charge)
- No consideration of plex/channel information
- All plexes share the same pairing relationships

### 3. How Channel Information is Used

#### In the Repository:
- Channel is loaded from parquet files (or set to 0 if missing)
- Channel is used for protein rollup grouping: `groupby(gdf, [:file_name, :channel, :decoy, :entrapment_group, :protein])`
- Channel is NOT used in EFDR calculations

#### In the Notebook:
- Channel (plex) is a core part of the pairing logic
- Each plex has independent original/entrapment score comparisons
- EFDR is calculated considering plex-specific pairs

### 4. Implications of the Difference

#### Scenario: Same peptide in multiple plexes
Consider a peptide "PEPTIDE" with charge 2 that appears in plexes 1, 2, and 3:

**Notebook approach:**
- Plex 1: (PEPTIDE, z=2) original score: 0.9, entrapment score: 0.7
- Plex 2: (PEPTIDE, z=2) original score: 0.6, entrapment score: 0.8
- Plex 3: (PEPTIDE, z=2) original score: 0.95, entrapment score: 0.5

Each plex independently determines if the entrapment "wins" or "loses".

**Repository approach:**
- Only one pairing exists for (PEPTIDE, z=2) across all plexes
- The complement_indices point to a single entrapment partner
- Cannot handle plex-specific score variations

### 5. Why This Matters

1. **Biological Relevance**: In multiplex experiments, the same peptide can have different abundances/scores across channels due to:
   - Different sample conditions
   - Technical variation
   - Biological variation between samples

2. **Statistical Impact**: The paired EFDR method specifically compares original vs entrapment scores. If these comparisons are not plex-specific:
   - We lose the ability to detect plex-specific false discoveries
   - The EFDR calculation becomes less accurate for multiplex data

3. **Current Repository Behavior**: 
   - For single-plex experiments: Works correctly
   - For multiplex experiments: May produce incorrect EFDR values

### 6. Recommended Fix

To align the repository with the notebook's behavior, we need to:

1. **Modify compute_pairing_vectors** to include channel/plex information:
```julia
# Instead of:
key = (lib_sequences[i], lib_charges[i])

# Use:
key = (lib_sequences[i], lib_charges[i], lib_channels[i])
```

2. **Update the pairing logic** to work within plex groups:
```julia
# Group results by plex first
plex_groups = groupby(results_df, :channel)

# Compute pairing vectors for each plex separately
for (plex_key, plex_df) in pairs(plex_groups)
    pairing_info = compute_pairing_vectors(library_df, plex_df)
    # Calculate EFDR for this plex
end
```

3. **Alternative approach**: Modify the EFDR calculation to consider plex-specific scores while maintaining global pairing structure.

### 7. Questions for Clarification

1. **Is the current behavior intentional?** Perhaps the repository is designed for a different use case where plex-specific pairing is not needed?

2. **What is the expected input data?** 
   - Single-plex experiments only?
   - Multiplex with combined analysis?
   - Multiplex with plex-specific analysis?

3. **Performance considerations**: The notebook's approach requires more memory (storing scores for each plex-pair combination). Is this acceptable?

## Conclusion

The repository's current implementation does not match the notebook's plex-specific pairing logic. This could lead to incorrect EFDR calculations for multiplex experiments. The fix requires incorporating channel/plex information into the pairing system, either by modifying the pairing vectors or by calculating EFDR separately for each plex.