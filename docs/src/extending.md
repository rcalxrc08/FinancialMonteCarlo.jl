# Extending the package

Models and payoffs are open for extension: implement a struct and one method, and all metrics,
simulation methods, positions and portfolios work with it.

## A new model

1. Define a struct subtyping the abstract type that describes it (see the hierarchy in
   [Models](@ref)): an `AbstractMonteCarloEngine` if it knows nothing about rates and dividends,
   otherwise one of the `AbstractMonteCarloProcess` subtypes, holding an `underlying` field.
2. Implement `FinancialMonteCarlo.simulate!(S, model, rfCurve, mcConfig, T)` (engines: without
   `rfCurve`) filling the preallocated `Nsim × (Nstep + 1)` matrix `S`.

Reusing an engine for the random part is usually enough. For instance, Black-Scholes parametrised
by the variance instead of the volatility:

```@example extending
using FinancialMonteCarlo

struct BlackScholesVarianceProcess{num <: Number, U <: FinancialMonteCarlo.AbstractUnderlying} <: FinancialMonteCarlo.ItoProcess{num}
    v::num # variance
    underlying::U
end

function FinancialMonteCarlo.simulate!(S, model::BlackScholesVarianceProcess, rfCurve::FinancialMonteCarlo.AbstractZeroRateCurve, mcConfig::FinancialMonteCarlo.AbstractMonteCarloConfiguration, T::Number)
    μ = rfCurve.r - FinancialMonteCarlo.dividend(model) # risk neutral drift
    FinancialMonteCarlo.simulate!(S, GeometricBrownianMotion(sqrt(model.v), μ, model.underlying.S0), mcConfig, T)
    return nothing
end

model = BlackScholesVarianceProcess(0.04, Underlying(100.0, 0.01))
mcConfig = MonteCarloConfiguration(10_000, 30)
pricer(model, ZeroRate(0.02), mcConfig, EuropeanOption(1.0, 100.0))
```

A model built on an engine inherits everything the engine supports: antithetic variates, Sobol
sequences, time dependent rates (`μ` above is a curve when `rfCurve` is one) and, with CUDA.jl
loaded, the GPU. Nothing else needs to be written.

A model that needs its own time stepping can still support all of them with a single
`simulate!` body, built on the hooks the built-in models use (see `src/models/paths.jl`):

- `step_drifts(μ, dt, Nstep)` and `drift_of_step(drifts, j)`: the drift of each step, a single
  number for a constant rate and a vector only for a curve;
- `path_sides(X, mcConfig)`: the paths to fill with the sign of their normal draws, all of them,
  or the two halves with antithetic variates;
- `draw_buffer(mcConfig, type, N)` and `draw_normals!(Z, mcConfig)`: buffers and normal draws on
  the device of the configuration (CPU, or GPU with `CudaMode`).

Sobol sequences draw path by path, so they need a separate method dispatching on
`FinancialMonteCarlo.SerialSobolMonteCarloConfig`.

## A new payoff

1. Define a struct subtyping the payoff family (`EuropeanPayoff`, `PathDependentPayoff`,
   `AmericanPayoff`, ...), with the maturity in a field named `T`.
2. Implement `FinancialMonteCarlo.payout`: for European and American payoffs it takes the value
   of the underlying at one time, for path dependent payoffs the whole path.

The generic engine of the family takes care of picking the right time steps, discounting and,
for early exercise, of the Longstaff-Schwartz regression. A capped call:

```@example extending
struct CappedCall{num <: Number} <: FinancialMonteCarlo.EuropeanPayoff{num}
    T::num
    K::num
    cap::num
end

FinancialMonteCarlo.payout(ST::Number, p::CappedCall) = min(max(ST - p.K, zero(ST)), p.cap)

capped = CappedCall(1.0, 100.0, 10.0)
pricer(model, ZeroRate(0.02), mcConfig, capped)
```

The new payoff takes part in the contract algebra like the built-in ones:

```@example extending
pricer(model, ZeroRate(0.02), mcConfig, 2.0 * capped - EuropeanOption(1.0, 110.0))
```

## Automatic differentiation friendliness

To keep dual numbers working:

- keep numeric fields parametric (`CappedCall` uses one shared parameter, the built-in structs one
  per field) so that they can hold dual numbers;
- allocate with the type returned by `FinancialMonteCarlo.predict_output_type_zero`, never `Float64`;
- wrap assertions and other non differentiable code in `ChainRulesCore.@ignore_derivatives`.
