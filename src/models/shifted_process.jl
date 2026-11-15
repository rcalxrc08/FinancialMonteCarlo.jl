"""
    ShiftedProcess(model, α)

Shifted (displaced) version of any single name `model`: the model started from ``(1 - α) S_0`` plus a
deterministic shift ``α S_0`` growing at the risk neutral rate,
```math
S_t = X_t + α S_0 e^{\\int_0^t (r_s - d_s)\\,ds}, \\qquad X_0 = (1 - α) S_0,
```
where ``X`` follows `model`. The forward is the one of `model`; with ``0 ≤ α < 1`` the price stays
above ``α`` times the forward, with ``α < 0`` it can go below zero. `ShiftedProcess(BlackScholesProcess(σ, underlying), α)`
is the displaced diffusion, and `ShiftedProcess(LogNormalMixture(η, λ, underlying), α)` the
[`ShiftedLogNormalMixture`](@ref).

The spot and the dividends are those of the underlying of `model`, held by the shifted process (a bump of the spot,
as in [`delta`](@ref), moves the shift too); the paths of `model` are simulated by its own `simulate!`, so every
sampling, mode and number type it supports works shifted too.

# Arguments
- `model`: a single name model, e.g. [`BlackScholesProcess`](@ref), [`BachelierProcess`](@ref), [`HestonProcess`](@ref).
- `α`: shift, as a fraction of the spot price, smaller than 1.

[`FinancialMonteCarlo.get_parameters`](@ref) returns the parameters of `model` followed by `α`.

# Examples
```julia
model = ShiftedProcess(BlackScholesProcess(0.2, Underlying(100.0, 0.01)), 0.3)
pricer(model, ZeroRate(0.02), MonteCarloConfiguration(10_000, 30), EuropeanOption(1.0, 100.0))
```
"""
struct ShiftedProcess{M <: ScalarMonteCarloProcess, num <: Number, abstrUnderlying <: AbstractUnderlying, numtype <: Number} <: ScalarMonteCarloProcess{numtype}
    model::M
    α::num
    underlying::abstrUnderlying
    function ShiftedProcess(model::M, α::num, underlying::abstrUnderlying) where {numM <: Number, M <: ScalarMonteCarloProcess{numM}, num <: Number, abstrUnderlying <: AbstractUnderlying}
        ChainRulesCore.@ignore_derivatives @assert α < 1 "the shift must be smaller than the spot price"
        zero_typed = ChainRulesCore.@ignore_derivatives zero(numM) + zero(num)
        return new{M, num, abstrUnderlying, typeof(zero_typed)}(model, α, underlying)
    end
end
ShiftedProcess(model::ScalarMonteCarloProcess, α::Number) = ShiftedProcess(model, α, model.underlying)

export ShiftedProcess;

function simulate!(X, mcProcess::ShiftedProcess, rfCurve::AbstractZeroRateCurve, mcBaseData::AbstractMonteCarloConfiguration, T::Number)
    ChainRulesCore.@ignore_derivatives @assert T > 0
    S0 = mcProcess.underlying.S0
    d = dividend(mcProcess)
    α = mcProcess.α
    unshifted = @set mcProcess.model.underlying = Underlying(S0 * (1 - α), d)
    simulate!(X, unshifted.model, rfCurve, mcBaseData, T)
    add_shift!(X, α * S0, rfCurve.r - d, mcBaseData, T)
    return nothing
end
