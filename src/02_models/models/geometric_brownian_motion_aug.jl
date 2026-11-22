"""
    GeometricBrownianMotionVec(σ, μ, x0)

Engine simulating a geometric Brownian motion with time dependent drift. Built by
[`GeometricBrownianMotion`](@ref) when the drift is a curve.

# Arguments
- `σ`: volatility, must be positive.
- `μ`: drift curve.
- `x0`: initial value.
"""
struct GeometricBrownianMotionVec{num <: Number, C <: AbstractCurve, num2 <: Number, numtype <: Number} <: AbstractGeometricBrownianMotion{numtype}
    σ::num
    μ::C
    x0::num2
    function GeometricBrownianMotionVec(σ::num, μ::C, x0::num2) where {num <: Number, C <: AbstractCurve, num2 <: Number}
        ChainRulesCore.@ignore_derivatives check_domain(σ > 0, σ, "Volatility must be positive")
        zero_typed = ChainRulesCore.@ignore_derivatives zero(num) + zero(number_type(μ)) + zero(num2)
        return new{num, C, num2, typeof(zero_typed)}(σ, μ, x0)
    end
end

function GeometricBrownianMotion(σ::num, μ::C, x0::num2) where {num <: Number, C <: AbstractCurve, num2 <: Number}
    return GeometricBrownianMotionVec(σ, μ, x0)
end

export GeometricBrownianMotionVec;
