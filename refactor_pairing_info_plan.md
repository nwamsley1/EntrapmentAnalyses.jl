# Plan to Refactor compute_pairing_vectors to Modify DataFrame Directly

## Current Situation

Currently, `compute_pairing_vectors` returns a named tuple `pairing_info` with:
- `is_original`: Bool vector
- `pair_indices`: Int vector  
- `entrap_labels`: Int vector
- `complement_indices`: Int vector (deprecated)
- `complement_scores`: Float32 vector

The API then uses these vectors like:
```julia
pairing_info = compute_pairing_vectors(library_df, results_no_decoys)
# Later: pairing_info.is_original[run_indices]
```

## Proposed Change

Have `compute_pairing_vectors` directly add columns to the DataFrame:
- `:is_original` - Boolean indicating if peptide is original (true) or entrapment (false)
- `:pair_id` - Integer ID linking original/entrapment pairs from library
- `:entrap_label` - Integer label (0=original, 1=entrapment) for EFDR calculation
- `:complement_indices` - Deprecated, kept for compatibility (always -1)
- `:complement_score` - Float32 plex-specific score of the paired peptide

## Benefits

1. **Simpler API**: No need to track separate `pairing_info` object
2. **More intuitive**: Columns are part of the data they describe
3. **Easier subsetting**: When grouping by file, columns automatically subset
4. **Already partially implemented**: Current code already adds some columns

## Implementation Steps

### Step 1: Modify compute_pairing_vectors Signature and Documentation

Change from:
```julia
"""
    compute_pairing_vectors(library_df::DataFrame, results_df::DataFrame; kwargs...)

Compute plex-aware pairing vectors with file and plex-specific complement scores.

# Returns
Named tuple with:
- `is_original`: Bool vector indicating if each result is an original sequence
- `pair_indices`: Int vector with pair IDs for each result
- `entrap_labels`: Int vector with entrapment group IDs
- `complement_indices`: Int vector (deprecated, always -1)
- `complement_scores`: Float32 vector with plex-specific complement scores
"""
function compute_pairing_vectors(library_df::DataFrame, results_df::DataFrame; kwargs...)
```

To:
```julia
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

# Returns
- Modified `results_df` with new columns added
"""
function compute_pairing_vectors!(library_df::DataFrame, results_df::DataFrame; kwargs...)
```

### Step 2: Update compute_pairing_vectors Implementation

Current implementation already adds columns via `add_plex_complement_scores!`:
- `:complement_score`
- `:is_original` 
- `:pair_id`

Need to also add:
- `:entrap_label` (computed from `:is_original`)
- `:complement_indices` (for compatibility, fill with -1)

```julia
function compute_pairing_vectors!(
    library_df::DataFrame,
    results_df::DataFrame;
    # ... same kwargs
)
    println("Computing plex-aware pairing information...")
    
    # ... existing implementation ...
    
    # Add entrap_label column
    println("Adding entrap_label column (0=original, 1=entrapment)...")
    results_df[!, :entrap_label] = [orig ? 0 : 1 for orig in results_df[!, :is_original]]
    
    # Add complement_indices for backward compatibility
    results_df[!, :complement_indices] = fill(-1, nrow(results_df))
    
    println("Added pairing columns: is_original, pair_id, entrap_label, complement_score, complement_indices")
    
    return results_df
end
```

### Step 3: Update API Usage in run_efdr_analysis with Documentation

#### Before:
```julia
# Step 7: Compute pairing vectors
pairing_info = compute_pairing_vectors(library_df, results_no_decoys; 
                                      show_progress=show_progress)
```

#### After:
```julia
# Step 7: Add pairing information to dataframe
# This adds columns for tracking original/entrapment pairs and their scores:
# - is_original: Boolean flag (true=original, false=entrapment)
# - pair_id: Links original/entrapment pairs from the library
# - entrap_label: Integer label for EFDR calculation (0=original, 1=entrapment)
# - complement_score: Plex-specific score of the paired peptide
# - complement_indices: Deprecated, kept for compatibility
compute_pairing_vectors!(library_df, results_no_decoys; 
                        show_progress=show_progress)
```

### Step 4: Update run_efdr_analysis Docstring

Add to the docstring:
```julia
"""
    run_efdr_analysis(parquet_files::Vector{String}, library_path::String; kwargs...)

Run empirical FDR analysis at the precursor level following the notebook implementation.

The function modifies the results DataFrame by adding pairing information columns:
- `is_original`: Boolean indicating if peptide is original (true) or entrapment (false)
- `pair_id`: Integer ID linking original/entrapment pairs
- `entrap_label`: Integer label (0=original, 1=entrapment) for EFDR
- `complement_score`: Plex-specific score of the paired peptide
- `complement_indices`: Deprecated column for backward compatibility

# Arguments
...
"""
```

### Step 5: Update API Usage Throughout

#### In run_efdr_analysis:
```julia
# Remove this line (~98):
entrap_labels = [orig ? 0 : 1 for orig in pairing_info.is_original]

# Change (~115):
run_entrap_labels = entrap_labels[run_indices]
# To:
run_entrap_labels = run_df[!, :entrap_label]

# Change (~135):
complement_scores = pairing_info.complement_scores[run_indices]
# To:
complement_scores = run_df[!, :complement_score]

# Change (~141):
pairing_info.is_original[run_indices]
# To:
run_df[!, :is_original]

# Change single file case (~180-181):
Float64.(pairing_info.complement_scores)
pairing_info.is_original
# To:
Float64.(results_no_decoys[!, :complement_score])
results_no_decoys[!, :is_original]
```

### Step 6: Update run_protein_efdr_analysis Similarly

Make analogous changes with same documentation approach.

### Step 7: Consider Backward Compatibility

Options:
1. Keep old `compute_pairing_vectors` as deprecated wrapper
2. Make this a breaking change with version bump
3. Provide migration guide in documentation

## Implementation Checklist

- [ ] Update compute_pairing_vectors docstring with column descriptions
- [ ] Modify compute_pairing_vectors to compute_pairing_vectors!
- [ ] Add entrap_label column creation with descriptive message
- [ ] Add complement_indices column (deprecated)
- [ ] Add detailed comment at Step 7 in api.jl explaining columns
- [ ] Update run_efdr_analysis docstring
- [ ] Update run_efdr_analysis to use columns instead of pairing_info
- [ ] Update run_protein_efdr_analysis docstring
- [ ] Update run_protein_efdr_analysis to use columns
- [ ] Remove all references to pairing_info variable
- [ ] Test with sample data
- [ ] Update CLAUDE.md documentation

## Example of Updated Code Flow

After implementation, the code will be cleaner:

```julia
# Add pairing information (modifies DataFrame in place)
compute_pairing_vectors!(library_df, results_no_decoys)

# Later, when processing by file:
for (key, run_df) in pairs(groupby(results_no_decoys, :file_name))
    # All pairing columns are automatically part of run_df
    paired_efdr = calculate_paired_efdr(
        run_df[!, :PredVal],
        run_df[!, :complement_score],  # Direct column access
        run_df[!, :is_original],       # No indexing needed
        run_df[!, :local_qvalue]
    )
end
```

This approach makes the data flow more transparent and easier to understand.