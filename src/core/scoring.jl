# Scoring utilities

"""
    monotonize!(values::AbstractVector{T}) where T<:AbstractFloat

Ensure FDR values are monotonically non-decreasing.
Works in-place by traversing backwards through sorted results.

# Algorithm
Starting from the end (worst scores), ensures each FDR is at least
as large as the FDR of better-scoring results.
"""
function monotonize!(values::AbstractVector{T}) where T<:AbstractFloat
    current_min = one(T)
    
    for i in length(values):-1:1
        if values[i] > current_min
            values[i] = current_min
        else
            current_min = values[i]
        end
    end
    
    return values
end

"""
    monotonize!(values::AbstractVector{Union{Missing, T}}) where T<:AbstractFloat

Handle missing values in monotonization.
"""
function monotonize!(values::AbstractVector{Union{Missing, T}}) where T<:AbstractFloat
    current_min = one(T)
    
    for i in length(values):-1:1
        if ismissing(values[i])
            continue
        end
        
        if values[i] > current_min
            values[i] = current_min
        else
            current_min = values[i]
        end
    end
    
    return values
end

"""
    monotonize(values::AbstractVector)

Non-mutating version of monotonize!
"""
function monotonize(values::AbstractVector)
    return monotonize!(copy(values))
end