# Abstract type for the geometric Brownian motion engines (constant and time dependent drift)
abstract type AbstractGeometricBrownianMotion{T <: Number} <: AbstractMonteCarloEngine{T} end
"""
    GeometricBrownianMotion(σ, μ, x0)

Engine simulating the geometric Brownian motion ``dX_t = μ X_t\\,dt + σ X_t\\,dW_t``, ``X_0 = x_0``.
If `μ` is a curve a [`GeometricBrownianMotionVec`](@ref) with time dependent drift is returned.

# Arguments
- `σ`: volatility, must be positive.
- `μ`: drift.
- `x0`: initial value.
"""
struct GeometricBrownianMotion{num <: Number, num1 <: Number, num2 <: Number, numtype <: Number} <: AbstractGeometricBrownianMotion{numtype}
    σ::num
    μ::num1
    x0::num2
    function GeometricBrownianMotion(σ::num, μ::num1, x0::num2) where {num <: Number, num1 <: Number, num2 <: Number}
        ChainRulesCore.@ignore_derivatives check_domain(σ > 0, σ, "Volatility must be positive")
        zero_typed = ChainRulesCore.@ignore_derivatives zero(num) + zero(num1) + zero(num2)
        return new{num, num1, num2, typeof(zero_typed)}(σ, μ, x0)
    end
end

export GeometricBrownianMotion;

function simulate!(X, mcProcess::AbstractGeometricBrownianMotion, mcBaseData::AbstractMonteCarloConfiguration, T::Number)
    ChainRulesCore.@ignore_derivatives check_domain(T > 0, T, "Time to Maturity must be positive")
    σ = mcProcess.σ
    μ = mcProcess.μ - σ^2 / 2
    simulate!(X, BrownianMotion(σ, μ), mcBaseData, T)
    exp_paths!(X, mcProcess.x0)

    nothing
end
