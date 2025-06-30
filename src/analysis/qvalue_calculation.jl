# Q-value calculation functions

"""
    validate_sort_order(values::AbstractVector, ascending::Bool=true)

Validate that a vector is properly sorted.
Returns true if sorted correctly, false otherwise.
"""
function validate_sort_order(values::AbstractVector, ascending::Bool=true)
    if ascending
        return issorted(values)
    else
        return issorted(values, rev=true)
    end
end

"""
    calculate_qvalues(scores::AbstractVector{T}, is_decoy::AbstractVector{Bool}) where T<:AbstractFloat

Calculate q-values using target-decoy approach.

# Arguments
- `scores`: Score values (higher is better)
- `is_decoy`: Boolean vector indicating decoy status

# Returns
- Vector of q-values
"""
function calculate_qvalues(
    scores::AbstractVector{T},
    is_decoy::AbstractVector{Bool}
) where T<:AbstractFloat
    n = length(scores)
    
    if length(is_decoy) != n
        error("Scores and is_decoy must have the same length")
    end
    
    # Create indices for sorting
    #indices = collect(1:n)
    
    # Sort by score (descending), targets before decoys
    #sort_order = sortperm(indices) do i
    #    # Primary: score (descending)
    #    # Secondary: target before decoy
    #    (-scores[i], is_decoy[i])
    #end
    sort_order = sortperm(collect(zip(scores, is_decoy)), by=x -> (-x[1], x[2]))
    # Validate sort order
    sorted_scores = scores[sort_order]
    if !validate_sort_order(sorted_scores, false)  # Should be descending
        @warn "Scores may not be properly sorted. This could affect q-value calculation."
    end
    
    # Calculate q-values
    qvals = zeros(T, n)
    Nτ, Nd = 0, 0
    
    for i in 1:n
        idx = sort_order[i]
        
        if is_decoy[idx]
            Nd += 1
        else
            Nτ += 1
        end
        
        if Nτ > 0
            qvals[idx] = Nd / Nτ
        else
            qvals[idx] = zero(T)
        end
    end
    
    # Monotonize - work on sorted order
    # Create sorted view of qvals
    sorted_qvals = qvals[sort_order]
    monotonize!(sorted_qvals)
    
    # Map back to original positions
    for (i, idx) in enumerate(sort_order)
        qvals[idx] = sorted_qvals[i]
    end
    
    # Validate monotonicity
    if !validate_sort_order(sorted_qvals, true)
        error("Q-values are not properly monotonic after monotonization")
    end
    
    return qvals
end

"""
    calculate_qvalues!(df::AbstractDataFrame; score_col=:PredVal, decoy_col=:decoy, sort_df=true)

Calculate local q-values and add to DataFrame.

# Arguments
- `df`: DataFrame to process
- `score_col`: Score column name
- `decoy_col`: Decoy column name  
- `sort_df`: Whether to sort the DataFrame (default: true)
"""
function calculate_qvalues!(
    df::AbstractDataFrame;
    score_col::Symbol = :PredVal,
    decoy_col::Symbol = :decoy,
    sort_df::Bool = true
)
    println("Calculating local q-values...")
    
    # Validate columns
    if !hasproperty(df, score_col)
        error("Score column $score_col not found")
    end
    if !hasproperty(df, decoy_col)
        error("Decoy column not found")
    end
    
    # Sort dataframe if requested
    if sort_df
        sort!(df, [score_col, decoy_col], rev=[true, false])
        
        # Validate sort
        if !validate_sort_order(df[!, score_col], false)
            @warn "DataFrame may not be properly sorted by score"
        end
    end
    
    # Calculate q-values
    qvals = calculate_qvalues(
        Float32.(df[!, score_col]),
        df[!, decoy_col]
    )
    
    # Add to dataframe
    df[!, :local_qvalue] = qvals
    
    # Print summary
    Nd = sum(df[!, decoy_col])
    Nτ = nrow(df) - Nd
    println("Local q-values: $(Nd) decoys, $(Nτ) targets")
    
    return nothing
end

"""
    calculate_global_qvalues!(df::AbstractDataFrame; score_col=:PredVal)

Calculate global q-values (best per precursor).
Adds :global_qvalue column to the dataframe.
"""
function calculate_global_qvalues!(
    df::AbstractDataFrame;
    score_col::Symbol = :PredVal
)
    println("Calculating global q-values...")
    
    # Group by precursor
    gdf = groupby(df, [:decoy, :stripped_seq, :z])
    
    # Get best scoring per group
    global_summary = combine(gdf) do group
        idx = argmax(group[:, score_col])
        return group[idx:idx, [:decoy, :stripped_seq, :z, score_col]]
    end
    
    # Sort for q-value calculation
    sort!(global_summary, [score_col, :decoy], rev=[true, false])
    
    # Validate sort
    if !validate_sort_order(global_summary[!, score_col], false)
        @warn "Global summary may not be properly sorted"
    end
    
    # Calculate q-values for the summary
    global_qvals = calculate_qvalues(
        Float32.(global_summary[!, score_col]),
        global_summary[!, :decoy]
    )
    
    global_summary[!, :global_qvalue] = global_qvals
    
    # Map back to original dataframe
    qval_dict = Dict(
        (row.decoy, row.stripped_seq, row.z) => row.global_qvalue
        for row in eachrow(global_summary)
    )
    
    df[!, :global_qvalue] = [qval_dict[(row.decoy, row.stripped_seq, row.z)] 
                              for row in eachrow(df)]
    
    # Print summary
    Nd = sum(global_summary[!, :decoy])
    Nτ = nrow(global_summary) - Nd
    println("Global q-values: $(Nd) decoys, $(Nτ) targets")
    
    return nothing
end

"""
    calculate_qvalues_per_file!(df::AbstractDataFrame; score_col=:PredVal, file_col=:file_name)

Calculate q-values separately for each file in the dataframe.
Adds :local_qvalue and :global_qvalue columns.

# Arguments
- `df`: DataFrame with PSM results
- `score_col`: Score column name (default: :PredVal)
- `file_col`: File identifier column (default: :file_name)
"""
function calculate_qvalues_per_file!(
    df::AbstractDataFrame;
    score_col::Symbol = :PredVal,
    file_col::Symbol = :file_name
)
    println("Calculating q-values per file...")
    
    # Initialize columns
    df[!, :local_qvalue] = zeros(Float32, nrow(df))
    df[!, :global_qvalue] = zeros(Float32, nrow(df))
    
    # Group by file
    if hasproperty(df, file_col)
        gdf = groupby(df, file_col)
        n_files = length(gdf)
        println("Processing $n_files files...")
        
        for (i, (key, file_df)) in enumerate(pairs(gdf))
            file_name = key[file_col]
            println("  File $i/$n_files: $file_name")
            
            # Calculate local q-values for this file
            calculate_qvalues!(file_df; score_col=score_col, sort_df=true)
            
            # Calculate global q-values for this file
            calculate_global_qvalues!(file_df; score_col=score_col)
            
            # Print file summary
            n_psms = nrow(file_df)
            n_decoys = sum(file_df.decoy)
            n_targets = n_psms - n_decoys
            println("    PSMs: $n_psms (targets: $n_targets, decoys: $n_decoys)")
        end
    else
        # Single file case
        println("No file column found, treating as single file...")
        calculate_qvalues!(df; score_col=score_col, sort_df=true)
        calculate_global_qvalues!(df; score_col=score_col)
    end
    
    return nothing
end