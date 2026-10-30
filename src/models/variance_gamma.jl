"""
    VarianceGammaProcess(σ, θ, κ, underlying)

Variance Gamma model: the log price is a Brownian motion with drift `θ` and volatility `σ`,
evaluated at a gamma subordinator with unit mean rate and variance rate `κ`. The drift is chosen to
make the discounted price a martingale.

# Arguments
- `σ`: volatility of the subordinated Brownian motion, must be positive.
- `θ`: drift of the subordinated Brownian motion.
- `κ`: variance rate of the subordinator, must be positive.
- `underlying`: an [`Underlying`](@ref) with spot price and dividends.

The parameters must satisfy ``1 - σ^2 κ / 2 - θ κ > 0`` for the forward to be finite.
"""
struct VarianceGammaProcess{num <: Number, num1 <: Number, num2 <: Number, abstrUnderlying <: AbstractUnderlying, numtype <: Number} <: InfiniteActivityProcess{numtype}
    σ::num
    θ::num1
    κ::num2
    underlying::abstrUnderlying
    function VarianceGammaProcess(σ::num, θ::num1, κ::num2, underlying::abstrUnderlying) where {num <: Number, num1 <: Number, num2 <: Number, abstrUnderlying <: AbstractUnderlying}
        ChainRulesCore.@ignore_derivatives @assert σ > 0 "volatility must be positive"
        ChainRulesCore.@ignore_derivatives @assert κ > 0 "κappa must be positive"
        ChainRulesCore.@ignore_derivatives @assert 1 - σ * σ * κ / 2 - θ * κ >= 0 "Parameters with unfeasible values"
        zero_typed = ChainRulesCore.@ignore_derivatives zero(num) + zero(num1) + zero(num2)
        return new{num, num1, num2, abstrUnderlying, typeof(zero_typed)}(σ, θ, κ, underlying)
    end
end

export VarianceGammaProcess;

compute_drift(::VarianceGammaProcess, σ::Number, θ::Number, κ::Number) = -1 / κ * log(1 - σ^2 * κ / 2 - θ * κ)
compute_drift(mcProcess::VarianceGammaProcess) = compute_drift(mcProcess, mcProcess.σ, mcProcess.θ, mcProcess.κ)
compute_subordinator_law(::VarianceGammaProcess, κ::Number, dt::Number) = Gamma(dt / κ, κ)
compute_subordinator_law(mcProcess::VarianceGammaProcess, dt::Number) = compute_subordinator_law(mcProcess, mcProcess.κ, dt)
