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
        ChainRulesCore.@ignore_derivatives @assert σ > 0 "Volatility must be positive"
        zero_typed = ChainRulesCore.@ignore_derivatives zero(num) + zero(num1)
        return new{num, num1, typeof(zero_typed)}(σ, μ)
    end
end

export BrownianMotion;

# function simulate_path!(X, mcProcess::BrownianMotion, mcBaseData::MonteCarloConfiguration, T::numb) where {numb <: Number}
#     Nstep = mcBaseData.Nstep
#     σ = mcProcess.σ
#     μ = mcProcess.μ
# ChainRulesCore.@ignore_derivatives @assert T > 0
#     dt = T / Nstep
#     mean_bm = μ * dt
#     stddev_bm = σ * sqrt(dt)
#     isDualZero = mean_bm * stddev_bm * 0
#     @views X[1] = isDualZero
#     @inbounds for j = 1:Nstep
#         @views X[j+1] = X[j] + mean_bm + stddev_bm * randn(host_rng(mcBaseData))
#     end
#     return
# end
