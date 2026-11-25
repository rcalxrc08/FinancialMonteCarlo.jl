"""
    BachelierProcess(σ)

Bachelier (normal) model: the price discounted with the carry is an arithmetic Brownian motion,
```math
S_t = e^{\\int_0^t (r_s - d_s)\\,ds} \\left(S_0 + σ W_t\\right),
```
so the discounted price with dividends, ``S_0 + σ W_t``, is a martingale. ``S_T`` is normal, with
mean the forward ``F = S_0 e^{\\int_0^T (r - d)}`` and standard deviation ``σ \\sqrt{T} e^{\\int_0^T (r - d)}``:
European prices are those of the Bachelier formula with that normal volatility (``σ`` itself when
``r = d``, as on a future). Prices can be negative.

# Arguments
- `σ`: normal volatility, in units of the price, must be positive.

# Examples
```julia
model = BachelierProcess(20.0)
pricer(model, Underlying(100.0, 0.01), ZeroRate(0.02), MonteCarloConfiguration(10_000, 30), EuropeanOption(1.0, 100.0))
```
"""
struct BachelierProcess{num <: Number} <: ItoProcess{num}
    σ::num
    function BachelierProcess(σ::num) where {num <: Number}
        check_volatility(σ)
        return new{num}(σ)
    end
end

export BachelierProcess;

function simulate!(S, mcProcess::BachelierProcess, und::AbstractUnderlying, rfCurve::AbstractZeroRateCurve, mcBaseData::AbstractMonteCarloConfiguration, T::Number)
    check_maturity(T)
    σ = mcProcess.σ
    simulate!(S, BrownianMotion(σ, zero(σ)), mcBaseData, T)
    S0 = und.S0
    affine_carry!(S, S0, carry_factors(rfCurve.r - dividend(und), mcBaseData, T))
    return nothing
end
