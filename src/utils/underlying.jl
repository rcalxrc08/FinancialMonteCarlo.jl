###########Underlying Types
#Abstract
abstract type AbstractUnderlying end

#Scalar dividend
struct UnderlyingScalar{num <: Number, num2 <: Number} <: AbstractUnderlying
    S0::num
    d::num2
    function UnderlyingScalar(S0::num_, d::num_2 = 0.0) where {num_ <: Number, num_2 <: Number}
        ChainRulesCore.@ignore_derivatives @assert S0 >= zero(num_) "Underlying starting value must be positive"
        return new{num_, num_2}(S0, d)
    end
end
#Curve dividend
struct UnderlyingVec{num <: Number, num2 <: Number, num3 <: Number} <: AbstractUnderlying
    S0::num
    d::CurveType{num2, num3}
    function UnderlyingVec(S0::num_, d::CurveType{num2, num3}) where {num_ <: Number, num2 <: Number, num3 <: Number}
        ChainRulesCore.@ignore_derivatives @assert S0 >= zero(num_) "Underlying starting value must be positive"
        return new{num_, num2, num3}(S0, d)
    end
end

#Common Constructors
"""
    Underlying(S0, d = 0.0)

Spot price `S0` and continuous dividend yield `d` of an asset, held by every single name model.
`d` is either a number or a time dependent curve (see [`FinancialMonteCarlo.Curve`](@ref) and
[`FinancialMonteCarlo.ImpliedCurve`](@ref)). `S0` must be non negative.

# Examples
```julia
Underlying(100.0)
Underlying(100.0, 0.01)
Underlying(100.0, FinancialMonteCarlo.Curve([0.01, 0.02], 1.0))
```
"""
function Underlying(S0::num_, d::CurveType{num2, num3}) where {num_ <: Number, num2 <: Number, num3 <: Number}
    return UnderlyingVec(S0, d)
end
function Underlying(S0::num_, d::num_2 = 0.0) where {num_ <: Number, num_2 <: Number}
    return UnderlyingScalar(S0, d)
end

export Underlying;
