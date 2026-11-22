###########Underlying Types
#Abstract
abstract type AbstractUnderlying end

#Scalar dividend
struct UnderlyingScalar{num <: Number, num2 <: Number} <: AbstractUnderlying
    S0::num
    d::num2
    function UnderlyingScalar(S0::num_, d::num_2 = zero(S0)) where {num_ <: Number, num_2 <: Number}
        check_spot(S0)
        return new{num_, num_2}(S0, d)
    end
end
#Curve dividend
struct UnderlyingVec{num <: Number, C <: Curve} <: AbstractUnderlying
    S0::num
    d::C
    function UnderlyingVec(S0::num_, d::C) where {num_ <: Number, C <: Curve}
        check_spot(S0)
        return new{num_, C}(S0, d)
    end
end

#Common Constructors
"""
    Underlying(S0, d = zero(S0))

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
function Underlying(S0::Number, d::Curve)
    return UnderlyingVec(S0, d)
end
function Underlying(S0::num_, d::num_2 = zero(S0)) where {num_ <: Number, num_2 <: Number}
    return UnderlyingScalar(S0, d)
end

export Underlying;
