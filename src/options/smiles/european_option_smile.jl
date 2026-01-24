struct EuropeanOptionSmile{T1, T2, T3 <: AbstractArray{T2}, T4 <: AbstractArray{Bool}} <: AbstractSmile{EuropeanOption}
    T::T1
    K::T3
    isCall::T4
    function EuropeanOptionSmile(T::T1, K::T3, isCall::T4) where {T1, T2, T3 <: (AbstractArray{T2}), T4 <: (AbstractArray{Bool})}
        validate_abstract_smile(T, K, isCall)
        new{T1, T2, T3, T4}(T, K, isCall)
    end
end