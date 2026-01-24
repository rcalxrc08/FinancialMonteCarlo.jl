abstract type AbstractSmile{T} end

function validate_abstract_smile(T::T1, K::T3, isCall::T4) where {T1, T2, T3 <: (AbstractArray{T2}), T4 <: (AbstractArray{Bool})}
    @assert length(K) == length(isCall) "K and isCall must have the same length"
    @assert minimum(K) > 0 "Strikes must be positive"
    @assert T > 0 "Maturity must be positive"
end

maturity(smile::AbstractSmile) = smile.T

function smile_to_options(smile::AbstractSmile{T}) where {T}
    return T.(smile.T, smile.K, smile.isCall)
end