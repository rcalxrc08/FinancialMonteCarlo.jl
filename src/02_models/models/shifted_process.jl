"""
    ShiftedProcess(model, α)

Shifted (displaced) version of any single name `model`: the model started from ``(1 - α) S_0`` plus a
deterministic shift ``α S_0`` growing at the risk neutral rate,
```math
S_t = X_t + α S_0 e^{\\int_0^t (r_s - d_s)\\,ds}, \\qquad X_0 = (1 - α) S_0,
```
where ``X`` follows `model`. The forward is the one of `model`; with ``0 ≤ α < 1`` the price stays
above ``α`` times the forward, with ``α < 0`` it can go below zero. `ShiftedProcess(BlackScholesProcess(σ), α)`
is the displaced diffusion, and `ShiftedProcess(LogNormalMixture(η, λ), α)` the shifted lognormal
mixture (it replaces `ShiftedLogNormalMixture`, removed).

The spot and the dividends are those of the underlying passed with the model (a bump of the spot, as in
[`delta`](@ref), moves the shift too); `model` is simulated from `Underlying((1 - α) S0, d)`; the paths of `model` are simulated by its own `simulate!`, so every
sampling, mode and number type it supports works shifted too.

# Arguments
- `model`: a single name model, e.g. [`BlackScholesProcess`](@ref), [`BachelierProcess`](@ref), [`HestonProcess`](@ref).
- `α`: shift, as a fraction of the spot price, smaller than 1.

[`FinancialMonteCarlo.get_parameters`](@ref) returns the parameters of `model` followed by `α`.

# Examples
```julia
model = ShiftedProcess(BlackScholesProcess(0.2), 0.3)
pricer(model, Underlying(100.0, 0.01), ZeroRate(0.02), MonteCarloConfiguration(10_000, 30), EuropeanOption(1.0, 100.0))
```
"""
struct ShiftedProcess{M <: ScalarMonteCarloProcess, num <: Number, numtype <: Number} <: ScalarMonteCarloProcess{numtype}
    model::M
    α::num
    function ShiftedProcess(model::M, α::num) where {numM <: Number, M <: ScalarMonteCarloProcess{numM}, num <: Number}
        check_shift(α)
        zero_typed = zero_of_types(numM, num)
        return new{M, num, typeof(zero_typed)}(model, α)
    end
end
check_shift(α) = check_domain(α < 1, α, "the shift must be smaller than the spot price")
ChainRulesCore.@non_differentiable check_shift(::Any)

export ShiftedProcess;

function simulate!(X, mcProcess::ShiftedProcess, und::AbstractUnderlying, rfCurve::AbstractZeroRateCurve, mcBaseData::AbstractMonteCarloConfiguration, T::Number)
    check_maturity(T)
    S0 = und.S0
    d = dividend(und)
    α = mcProcess.α
    simulate!(X, mcProcess.model, Underlying(S0 * (1 - α), d), rfCurve, mcBaseData, T)
    add_shift!(X, α * S0, rfCurve.r - d, mcBaseData, T)
    return nothing
end
