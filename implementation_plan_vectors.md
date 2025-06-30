# Implementation Plan Update: Vector-Based API Design

## Core Design Change

Instead of passing DataFrames and column symbols to functions, we'll pass abstract vectors directly. This provides:
- Type stability and better performance
- Clearer function signatures
- Easier testing with mock data
- Decoupling from specific column names

## Updated Function Signatures

### Pairing System
```julia
# Core function that works with vectors
function compute_pairing_vectors(
    lib_sequences::AbstractVector{String},
    lib_charges::AbstractVector{Int},
    lib_entrap_groups::AbstractVector{Int},
    lib_pair_ids::AbstractVector{Int},
    results_sequences::AbstractVector{String},
    results_charges::AbstractVector{Int};
    show_progress::Bool = true
)
    # Returns named tuple of vectors
end

# Convenience wrapper for DataFrames
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
    # Extract vectors and call core function
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
```

### EFDR Methods
```julia
# Core EFDR calculation with vectors only
function calculate_paired_efdr(
    scores::AbstractVector{T},
    complement_scores::AbstractVector{T},
    is_original::AbstractVector{Bool},
    qvals::AbstractVector{T};
    r::T = one(T),
    show_progress::Bool = true
) where T<:AbstractFloat
    # Direct calculation on vectors
end

# DataFrame wrapper
function calculate_paired_efdr(
    df::DataFrame;
    score_col::Symbol = :PredVal,
    qval_col::Symbol = :local_qvalue,
    pairing_info::NamedTuple,
    r::Float64 = 1.0,
    show_progress::Bool = true
)
    # Extract complement scores
    complement_scores = get_complement_scores(
        df[!, score_col], 
        pairing_info.complement_indices
    )
    
    # Call vector-based function
    return calculate_paired_efdr(
        df[!, score_col],
        complement_scores,
        pairing_info.is_original,
        df[!, qval_col];
        r = r,
        show_progress = show_progress
    )
end
```

### Q-value Calculations
```julia
# Core function with vectors
function calculate_qvalues(
    scores::AbstractVector{T},
    is_decoy::AbstractVector{Bool}
) where T<:AbstractFloat
    # Returns q-values vector
end

# DataFrame wrapper
function calculate_qvalues!(
    df::DataFrame;
    score_col::Symbol = :PredVal,
    decoy_col::Symbol = :decoy
)
    qvals = calculate_qvalues(
        df[!, score_col],
        df[!, decoy_col]
    )
    df[!, :local_qvalue] = qvals
    return nothing
end
```

## Benefits of This Approach

1. **Performance**: Julia's compiler can optimize vector operations better
2. **Testing**: Easy to test with simple vector inputs
3. **Flexibility**: Column names only matter at the API boundary
4. **Type Safety**: Function signatures clearly show expected types
5. **Reusability**: Core functions can be used outside DataFrame context

## Implementation Strategy

1. Implement all core algorithms as vector-based functions
2. Add thin DataFrame wrappers that handle column extraction
3. Keep column name mappings in one place (API functions)
4. Use abstract types to allow different vector implementations

This design follows Julia best practices and makes the code more maintainable and performant.