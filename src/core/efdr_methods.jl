# EFDR calculation methods

abstract type EFDRMethod end

"""
    calculate_combined_efdr(scores, entrap_labels, qvals; r=1.0, show_progress=true)

Calculate combined empirical FDR using vector inputs.

# Arguments
- `scores`: Score values
- `entrap_labels`: Entrapment group labels (0 = original, >0 = entrapment)
- `qvals`: Q-values for sorting
- `r`: Ratio parameter (default: 1.0)
- `show_progress`: Show progress bar (default: true)

# Returns
- Vector of empirical FDR values
"""
function calculate_combined_efdr(
    scores::AbstractVector{T},
    entrap_labels::AbstractVector{<:Integer},
    qvals::AbstractVector{T};
    r::T = one(T),
    show_progress::Bool = true
) where T<:AbstractFloat
    n = length(scores)
    
    # Validate inputs
    if length(entrap_labels) != n || length(qvals) != n
        error("All input vectors must have the same length")
    end
    
    efdr = zeros(T, n)
    
    # Sort by q-value (ascending), then score (descending)
    sort_order = sortperm(collect(zip(qvals, -scores)))
    
    # Validate sort order
    sorted_qvals = qvals[sort_order]
    if !issorted(sorted_qvals)
        @warn "Q-values are not properly sorted. This may affect EFDR calculation accuracy."
    end
    
    Nτ, Nε = 0, 0
    
    pb = show_progress ? ProgressBar(1:n) : (1:n)
    if show_progress
        set_description(pb, "Calculating Combined EFDR")
    end
    
    for i in pb
        idx = sort_order[i]
        
        if entrap_labels[idx] == 0
            Nτ += 1
        else
            Nε += 1
        end
        
        if Nτ + Nε > 0
            efdr[idx] = min(one(T), (Nε * (1 + 1/r)) / (Nτ + Nε))
        end
    end
    
    return efdr
end

"""
    calculate_paired_efdr(scores, complement_scores, is_original, qvals; r=1.0, show_progress=true)

Calculate paired empirical FDR using vector inputs with O(n²) complexity.

# Arguments
- `scores`: Score values for all sequences
- `complement_scores`: Scores of paired sequences (-1 if no pair)
- `is_original`: Boolean vector indicating original (true) vs entrapment (false)
- `qvals`: Q-values for sorting
- `r`: Ratio parameter (default: 1.0)
- `show_progress`: Show progress bar (default: true)

# Returns
- Vector of empirical FDR values
"""
function calculate_paired_efdr(
    scores::AbstractVector{T},
    complement_scores::AbstractVector{T},
    is_original::AbstractVector{Bool},
    qvals::AbstractVector{T};
    r::T = one(T),
    show_progress::Bool = true
) where T<:AbstractFloat
    n = length(scores)
    
    # Validate inputs
    if length(complement_scores) != n || length(is_original) != n || length(qvals) != n
        error("All input vectors must have the same length")
    end
    
    efdr = zeros(T, n)
    
    # Sort by q-value (ascending), then score (descending)
    sort_order = sortperm(collect(zip(qvals, -scores)))
    
    # Validate sort order
    sorted_qvals = qvals[sort_order]
    if !issorted(sorted_qvals)
        @warn "Q-values are not properly sorted. This may affect EFDR calculation accuracy."
    end
    
    # Pre-compute total operations for accurate progress
    total_ops = sum(1:n)  # n*(n+1)/2
    completed_ops = 0
    
    # Pre-compute score relationships
    println("Pre-computing score relationships...")
    score_rels = Vector{Symbol}(undef, n)
    
    for i in 1:n
        if is_original[i]
            score_rels[i] = :original
        else
            e_score = scores[i]
            o_score = complement_scores[i]
            
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
        s = scores[sort_order[i]]
        
        for j in 1:i
            idx = sort_order[j]
            
            if is_original[idx]
                Nτ += 1
            else
                Nε += 1
                
                if score_rels[idx] == :entrap_wins
                    o_score = complement_scores[idx]
                    if o_score >= s
                        Nετs += 1
                    end
                elseif score_rels[idx] != :unpaired && scores[idx] >= s
                    o_score = complement_scores[idx]
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
            efdr[sort_order[i]] = min(one(T), (Nε + Nεsτ + 2*Nετs) / (Nτ + Nε))
        end
    end
    
    if show_progress
        finish!(pb)
    end
    
    return efdr
end

# Struct-based interface for compatibility with Pioneer.jl style

struct CombinedEFDR{T<:Real} <: EFDRMethod
    scores::Vector{T}
    entrap_labels::Vector{Int}
    qvals::Vector{T}
    r::T
end

struct PairedEFDR{T<:Real} <: EFDRMethod
    scores::Vector{T}
    complement_scores::Vector{T}
    is_original::Vector{Bool}
    qvals::Vector{T}
    r::T
end

function calculate_efdr(method::CombinedEFDR; show_progress::Bool=true)
    return calculate_combined_efdr(
        method.scores,
        method.entrap_labels,
        method.qvals;
        r = method.r,
        show_progress = show_progress
    )
end

function calculate_efdr(method::PairedEFDR; show_progress::Bool=true)
    # Create entrap_labels from is_original for compatibility
    entrap_labels = [orig ? 0 : 1 for orig in method.is_original]
    
    return calculate_paired_efdr(
        method.scores,
        method.complement_scores,
        method.is_original,
        method.qvals;
        r = method.r,
        show_progress = show_progress
    )
end

# DataFrame convenience functions

"""
    calculate_paired_efdr(df::DataFrame; kwargs...)

DataFrame wrapper for paired EFDR calculation.

# Keyword Arguments
- `score_col`: Score column name (default: :PredVal)
- `qval_col`: Q-value column name (default: :local_qvalue)
- `pairing_info`: Named tuple from compute_pairing_vectors
- `r`: Ratio parameter (default: 1.0)
- `show_progress`: Show progress bar (default: true)
"""
function calculate_paired_efdr(
    df::DataFrame;
    score_col::Symbol = :PredVal,
    qval_col::Symbol = :local_qvalue,
    pairing_info::NamedTuple,
    r::Float64 = 1.0,
    show_progress::Bool = true
)
    # Get complement scores
    complement_scores = get_complement_scores(
        df[!, score_col],
        pairing_info.complement_indices
    )
    
    return calculate_paired_efdr(
        Float64.(df[!, score_col]),
        Float64.(complement_scores),
        pairing_info.is_original,
        Float64.(df[!, qval_col]);
        r = r,
        show_progress = show_progress
    )
end