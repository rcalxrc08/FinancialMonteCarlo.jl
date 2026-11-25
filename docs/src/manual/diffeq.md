```@raw html
<div class="fc-page-head">
<p class="fc-kicker"><span>Manual 11</span> Any SDE as a model</p>
```

# StochasticDiffEq.jl

When `StochasticDiffEq` is loaded, any SDE (or jump diffusion, with JumpProcesses) defined with the SciML problem types can
be used as a model through `FinancialMonteCarlo.MonteCarloDiffEqModel`. DifferentialEquations 7 loads StochasticDiffEq
itself; DifferentialEquations 8 contains only the ODE solvers, so with it load StochasticDiffEq (and JumpProcesses for
jumps) explicitly:

```@raw html
<div class="fc-meta"><span>StochasticDiffEq.jl</span><span>JumpProcesses.jl</span></div>
<nav class="fc-toc">
<a href="#Black-Scholes-as-an-SDE"><b>01</b>Black-Scholes as an SDE</a>
<a href="#A-jump-diffusion-in-log-price"><b>02</b>A jump diffusion in log price</a>
</nav>
</div>
```

```julia
FinancialMonteCarlo.MonteCarloDiffEqModel(ensemble_problem)
FinancialMonteCarlo.MonteCarloDiffEqModel(ensemble_problem, final_transform)
```

`ensemble_problem` is an `EnsembleProblem` (SciMLBase, reexported by StochasticDiffEq); the optional `final_transform`
is applied to the solution, e.g. to go from log prices to prices. The problem is solved with
`SOSRI()` and a fixed step equal to the time step of the [`MonteCarloConfiguration`](@ref), one
trajectory per path (`EnsembleThreads()`), and then priced like any other model. The zero rate
curve must be a scalar [`ZeroRate`](@ref), and the drift of the SDE must already be the risk
neutral one.

## Black-Scholes as an SDE

```julia
using FinancialMonteCarlo, StochasticDiffEq

S0, K, r, T, d, σ = 100.0, 100.0, 0.02, 1.0, 0.01, 0.2
drift(u, p, t) = (r - d) * u
diffusion(u, p, t) = σ * u
problem = EnsembleProblem(SDEProblem(drift, diffusion, S0, (0.0, T)))
model = FinancialMonteCarlo.MonteCarloDiffEqModel(problem)

pricer(model, Underlying(S0, d), ZeroRate(r), MonteCarloConfiguration(10_000, 30), EuropeanOption(T, K))
```

## A jump diffusion in log price

The Kou model, written as an SDE for the log price with a `ConstantRateJump`:

```julia
using JumpProcesses, Random
mcConfig = MonteCarloConfiguration(10_000, 30)
rng = FinancialMonteCarlo.host_rng(mcConfig)
p, λ, λ₊, λ₋ = 0.3, 5.0, 30.0, 20.0

drift(u, p_, t) = r - d - σ^2 / 2 - λ * (p / (λ₊ - 1) - (1 - p) / (λ₋ + 1))
diffusion(u, p_, t) = σ
rate(u, p_, t) = λ * T
jump!(integrator) = (integrator.u += rand(rng) < p ? randexp(rng) / λ₊ : -randexp(rng) / λ₋)

jump_problem = JumpProblem(SDEProblem(drift, diffusion, 0.0, (0.0, T)), Direct(), ConstantRateJump(rate, jump!))
model = FinancialMonteCarlo.MonteCarloDiffEqModel(EnsembleProblem(jump_problem), x -> S0 * exp(x))

pricer(model, Underlying(S0, d), ZeroRate(r), mcConfig, EuropeanOption(T, K))
```

The jump sizes draw from the generator of the configuration, which `set_seed!` resets; the jump times and the Brownian
increments come from the solver's own generator.
