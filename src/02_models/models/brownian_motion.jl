"""
    BrownianMotion(σ, μ)

Engine simulating the arithmetic Brownian motion ``X_t = μ t + σ W_t``, ``X_0 = 0``.
If `μ` is a curve (see [`FinancialMonteCarlo.Curve`](@ref)) a [`FinancialMonteCarlo.BrownianMotionVec`](@ref)
with time dependent drift is returned.

# Arguments
- `σ`: volatility, must be positive.
- `μ`: drift.
"""
struct BrownianMotion{num <: Number, num1 <: Number, numtype <: Number} <: AbstractMonteCarloEngine{numtype}
    σ::num
    μ::num1
    function BrownianMotion(σ::num, μ::num1) where {num <: Number, num1 <: Number}
        check_volatility(σ)
        zero_typed = zero_of_types(num, num1)
        return new{num, num1, typeof(zero_typed)}(σ, μ)
    end
end

export BrownianMotion;
