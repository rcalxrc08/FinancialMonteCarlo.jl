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
        ChainRulesCore.@ignore_derivatives @assert σ > 0 "Volatility must be positive"
        return new{num, abstrUnderlying}(σ, underlying)
    end
end

export BlackScholesProcess;

function simulate!(S, mcProcess::BlackScholesProcess, rfCurve::AbstractZeroRateCurve, mcBaseData::AbstractMonteCarloConfiguration, T::Number)
    ChainRulesCore.@ignore_derivatives @assert T > 0.0
    r = rfCurve.r
    d = dividend(mcProcess)
    σ = mcProcess.σ
    μ = r - d

    simulate!(S, GeometricBrownianMotion(σ, μ, mcProcess.underlying.S0), mcBaseData, T)

    nothing
end

# function simulate_path!(S, mcProcess::BlackScholesProcess, rfCurve::AbstractZeroRateCurve, mcBaseData::AbstractMonteCarloConfiguration, T::Number)
# ChainRulesCore.@ignore_derivatives @assert T > 0.0
#     r = rfCurve.r
#     d = dividend(mcProcess)
#     σ = mcProcess.σ
#     μ = r - d

#     simulate_path!(S, GeometricBrownianMotion(σ, μ, mcProcess.underlying.S0), mcBaseData, T)

#     nothing
# end
