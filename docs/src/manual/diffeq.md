# DifferentialEquations.jl

When `DifferentialEquations` is loaded, any SDE (or jump diffusion) defined with it can be used as
a model through `FinancialMonteCarlo.MonteCarloDiffEqModel`:

```julia
FinancialMonteCarlo.MonteCarloDiffEqModel(ensemble_problem, underlying)
FinancialMonteCarlo.MonteCarloDiffEqModel(ensemble_problem, final_transform, underlying)
```

`ensemble_problem` is a `DifferentialEquations.EnsembleProblem`; the optional `final_transform`
is applied to the solution, e.g. to go from log prices to prices. The problem is solved with
`SOSRI()` and a fixed step equal to the time step of the [`MonteCarloConfiguration`](@ref), one
trajectory per path (`EnsembleThreads()`), and then priced like any other model. The zero rate
curve must be a scalar [`ZeroRate`](@ref), and the drift of the SDE must already be the risk
neutral one.

## Black-Scholes as an SDE

```julia
using FinancialMonteCarlo, DifferentialEquations

S0, K, r, T, d, σ = 100.0, 100.0, 0.02, 1.0, 0.01, 0.2
drift(u, p, t) = (r - d) * u
diffusion(u, p, t) = σ * u
problem = EnsembleProblem(SDEProblem(drift, diffusion, S0, (0.0, T)))
model = FinancialMonteCarlo.MonteCarloDiffEqModel(problem, Underlying(S0, d))

pricer(model, ZeroRate(r), MonteCarloConfiguration(10_000, 30), EuropeanOption(T, K))
```

## A jump diffusion in log price

The Kou model, written as an SDE for the log price with a `ConstantRateJump`:

```julia
using Random
mcConfig = MonteCarloConfiguration(10_000, 30)
rng = mcConfig.parallelMode.rng
p, λ, λ₊, λ₋ = 0.3, 5.0, 30.0, 20.0

drift(u, p_, t) = r - d - σ^2 / 2 - λ * (p / (λ₊ - 1) - (1 - p) / (λ₋ + 1))
diffusion(u, p_, t) = σ
rate(u, p_, t) = λ * T
jump!(integrator) = (integrator.u += rand(rng) < p ? randexp(rng) / λ₊ : -randexp(rng) / λ₋)

jump_problem = JumpProblem(SDEProblem(drift, diffusion, 0.0, (0.0, T)), Direct(), ConstantRateJump(rate, jump!), rng = rng)
model = FinancialMonteCarlo.MonteCarloDiffEqModel(EnsembleProblem(jump_problem), x -> S0 * exp(x), Underlying(S0, d))

pricer(model, ZeroRate(r), mcConfig, EuropeanOption(T, K))
```

Using the generator of the configuration keeps the simulation reproducible.
