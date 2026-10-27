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

function simulate!(X, mcProcess::VarianceGammaProcess, rfCurve::AbstractZeroRateCurve, mcBaseData::AbstractMonteCarloConfiguration, T::Number)
    r = rfCurve.r
    S0 = mcProcess.underlying.S0
    d = dividend(mcProcess)
    Nstep = mcBaseData.Nstep
    σ = mcProcess.σ
    θ = mcProcess.θ
    κ = mcProcess.κ
    ChainRulesCore.@ignore_derivatives @assert T > 0
    dt = T / Nstep
    ψ = -1 / κ * log(1 - σ^2 * κ / 2 - θ * κ)
    drift = r - d - ψ

    gammaRandomVariable = Gamma(dt / κ, κ)

    simulate!(X, SubordinatedBrownianMotion(σ, drift, θ, gammaRandomVariable), mcBaseData, T)

    @. X = S0 * exp(X)

    nothing
end
