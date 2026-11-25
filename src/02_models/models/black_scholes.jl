"""
    BlackScholesProcess(σ)

Black-Scholes model: under the risk neutral measure the underlying follows a geometric Brownian
motion ``dS_t = (r_t - d_t) S_t\\,dt + σ S_t\\,dW_t``.

# Arguments
- `σ`: volatility, must be positive.

# Examples
```julia
model = BlackScholesProcess(0.2)
pricer(model, Underlying(100.0, 0.01), ZeroRate(0.02), MonteCarloConfiguration(10_000, 30), EuropeanOption(1.0, 100.0))
```
"""
struct BlackScholesProcess{num <: Number} <: ItoProcess{num}
    σ::num
    function BlackScholesProcess(σ::num) where {num <: Number}
        check_volatility(σ)
        return new{num}(σ)
    end
end

export BlackScholesProcess;

function simulate!(S, mcProcess::BlackScholesProcess, und::AbstractUnderlying, rfCurve::AbstractZeroRateCurve, mcBaseData::AbstractMonteCarloConfiguration, T::Number)
    check_maturity(T)
    r = rfCurve.r
    d = dividend(und)
    σ = mcProcess.σ
    μ = r - d

    simulate!(S, GeometricBrownianMotion(σ, μ, und.S0), mcBaseData, T)

    nothing
end
