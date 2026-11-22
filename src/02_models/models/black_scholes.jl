"""
    BlackScholesProcess(σ, underlying)

Black-Scholes model: under the risk neutral measure the underlying follows a geometric Brownian
motion ``dS_t = (r_t - d_t) S_t\\,dt + σ S_t\\,dW_t``.

# Arguments
- `σ`: volatility, must be positive.
- `underlying`: an [`Underlying`](@ref) with spot price and dividends.

# Examples
```julia
model = BlackScholesProcess(0.2, Underlying(100.0, 0.01))
pricer(model, ZeroRate(0.02), MonteCarloConfiguration(10_000, 30), EuropeanOption(1.0, 100.0))
```
"""
struct BlackScholesProcess{num <: Number, abstrUnderlying <: AbstractUnderlying} <: ItoProcess{num}
    σ::num
    underlying::abstrUnderlying
    function BlackScholesProcess(σ::num, underlying::abstrUnderlying) where {num <: Number, abstrUnderlying <: AbstractUnderlying}
        ChainRulesCore.@ignore_derivatives check_domain(σ > 0, σ, "Volatility must be positive")
        return new{num, abstrUnderlying}(σ, underlying)
    end
end

export BlackScholesProcess;

function simulate!(S, mcProcess::BlackScholesProcess, rfCurve::AbstractZeroRateCurve, mcBaseData::AbstractMonteCarloConfiguration, T::Number)
    ChainRulesCore.@ignore_derivatives check_domain(T > 0, T, "Time to Maturity must be positive")
    r = rfCurve.r
    d = dividend(mcProcess)
    σ = mcProcess.σ
    μ = r - d

    simulate!(S, GeometricBrownianMotion(σ, μ, mcProcess.underlying.S0), mcBaseData, T)

    nothing
end
