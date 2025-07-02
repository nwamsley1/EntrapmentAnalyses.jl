# Protein-level analysis functions

"""
    add_entrapment_group_column!(results_df::DataFrame, target_seqs::Set{String}, entrap_seqs::Set{String})

Add entrapment_group column to results DataFrame based on sequence membership.
"""
function add_entrapment_group_column!(
    results_df::DataFrame,
    target_seqs::Set{String},
    entrap_seqs::Set{String}
)
    # Type-stable vector operation
    stripped_seqs = results_df.stripped_seq::AbstractVector{String}
    entrapment_group = Vector{Bool}(undef, length(stripped_seqs))
    
    for i in eachindex(stripped_seqs)
        entrapment_group[i] = stripped_seqs[i] ∈ entrap_seqs
    end
    
    results_df[!, :entrapment_group] = entrapment_group
    return nothing
end

"""
    calculate_protein_qvalues_per_run!(protein_df::DataFrame)

Calculate protein-level q-values separately for each run.
"""
function calculate_protein_qvalues_per_run!(protein_df::DataFrame)
    protein_df[!, :Protein_Qvalue] = zeros(Float32, nrow(protein_df))
    
    for (run, results) in pairs(groupby(protein_df, :file_name))
        # Extract typed vectors
        decoy_vec = results.decoy::AbstractVector{Bool}
        n = length(decoy_vec)
        qvals = zeros(Float32, n)
        
        Nτ, Nd = 0, 0
        
        for i in 1:n
            if decoy_vec[i]
                Nd += 1
            else
                Nτ += 1
            end
            qvals[i] = Nd / Nτ
        end
        
        println("Run $(run.file_name): N_d $Nd N_t $Nτ")
        println("qvals $qvals")
        # Monotonize and assign back
        monotonize!(qvals)
        println("post monotonize! qvals $qvals")
        results[!, :Protein_Qvalue] = qvals
    end
    
    return nothing
end

"""
    prepare_protein_analysis(results_df::DataFrame, library_df::DataFrame; score_col=:PredVal)

Prepare data for protein-level analysis following the notebook approach.
"""
function prepare_protein_analysis(
    results_df::DataFrame,
    library_df::DataFrame;
    score_col::Symbol = :PredVal
)
    println("Preparing protein-level analysis...")
    
    # Step 1: Create target/entrapment sets from library
    entrap_groups = library_df.EntrapmentGroupId::AbstractVector{Int}
    peptide_seqs = library_df.PeptideSequence::AbstractVector{String}
    
    target_seqs = Set{String}()
    entrap_seqs = Set{String}()
    
    for i in eachindex(entrap_groups)
        if entrap_groups[i] == 0
            push!(target_seqs, peptide_seqs[i])
        else
            push!(entrap_seqs, peptide_seqs[i])
        end
    end
    
    # Step 2: Add entrapment_group column
    add_entrapment_group_column!(results_df, target_seqs, entrap_seqs)
    
    # Step 3: Sort by score (descending), targets above decoys for ties
    sort!(results_df, [score_col, :decoy], rev=[true, false])
    
    # Validate sort
    score_vec = results_df[!, score_col]::AbstractVector{Float32}
    if !issorted(score_vec, rev=true)
        @warn "Results may not be properly sorted by score"
    end
    
    # Step 4: Group and get best PSM per protein group
    gresults_df = groupby(results_df, [:file_name, :channel, :decoy, :entrapment_group, :protein])
    
    protein_df = combine(gresults_df) do group
        score_vec = group[!, score_col]::AbstractVector{Float32}
        idx = argmax(score_vec)
        return group[idx:idx, :]
    end
    
    # Step 5: Sort again for q-value calculation
    sort!(protein_df, [score_col, :decoy], rev=[true, false])
    
    # Step 6: Calculate protein q-values per run
    calculate_protein_qvalues_per_run!(protein_df)
    
    return protein_df
end


"""
    rollup_to_protein_groups(scores::AbstractVector{T}, proteins::AbstractVector{String},
                            file_names::AbstractVector{String}, channels::AbstractVector{<:Integer},
                            is_decoy::AbstractVector{Bool}, is_entrapment::AbstractVector{Bool}) where T

Vector-based protein rollup function.
"""
function rollup_to_protein_groups(
    scores::AbstractVector{T},
    proteins::AbstractVector{String},
    file_names::AbstractVector{String},
    channels::AbstractVector{<:Integer},
    is_decoy::AbstractVector{Bool},
    is_entrapment::AbstractVector{Bool}
) where T<:Real
    n = length(scores)
    
    # Validate inputs
    if length(proteins) != n || length(file_names) != n || length(channels) != n || 
       length(is_decoy) != n || length(is_entrapment) != n
        error("All input vectors must have the same length")
    end
    
    # Create unique groups
    groups = Dict{Tuple{String, Int, Bool, Bool, String}, Tuple{T, Int}}()
    
    for i in 1:n
        key = (file_names[i], channels[i], is_decoy[i], is_entrapment[i], proteins[i])
        
        if !haskey(groups, key) || scores[i] > groups[key][1]
            groups[key] = (scores[i], i)
        end
    end
    
    # Extract results
    n_groups = length(groups)
    protein_scores = Vector{T}(undef, n_groups)
    protein_files = Vector{String}(undef, n_groups)
    protein_channels = Vector{Int}(undef, n_groups)
    protein_is_decoy = Vector{Bool}(undef, n_groups)
    protein_is_entrapment = Vector{Bool}(undef, n_groups)
    protein_names = Vector{String}(undef, n_groups)
    original_indices = Vector{Int}(undef, n_groups)
    
    for (i, (key, (score, idx))) in enumerate(groups)
        protein_files[i] = key[1]
        protein_channels[i] = key[2]
        protein_is_decoy[i] = key[3]
        protein_is_entrapment[i] = key[4]
        protein_names[i] = key[5]
        protein_scores[i] = score
        original_indices[i] = idx
    end
    
    return (
        protein_scores = protein_scores,
        protein_names = protein_names,
        protein_files = protein_files,
        protein_channels = protein_channels,
        protein_is_decoy = protein_is_decoy,
        protein_is_entrapment = protein_is_entrapment,
        original_indices = original_indices
    )
end

"""
    rollup_to_protein_groups(df::DataFrame; score_col=:PredVal)

DataFrame wrapper for protein rollup.
"""
function rollup_to_protein_groups(df::DataFrame; score_col::Symbol = :PredVal)
    # Extract typed vectors
    scores = df[!, score_col]::AbstractVector{Float32}
    proteins = df.protein::AbstractVector{String}
    file_names = df.file_name::AbstractVector{String}
    channels = df.channel::AbstractVector{UInt8}
    is_decoy = df.decoy::AbstractVector{Bool}
    is_entrapment = df.entrapment_group::AbstractVector{Bool}
    
    return rollup_to_protein_groups(
        scores, proteins, file_names, Int.(channels), is_decoy, is_entrapment
    )
end