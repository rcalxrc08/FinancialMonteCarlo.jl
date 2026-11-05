```@raw html
<img class="pkg-hero-logo" src="assets/logo.png" alt="FinancialMonteCarlo.jl logo" width="96" height="96">
```

# FinancialMonteCarlo.jl

```@raw html
<p class="pkg-roles"><b>Monte Carlo pricing</b> · <b>Lévy processes</b> · <b>Automatic differentiation</b> · <b>GPU</b></p>
<div class="pkg-intro">
<p>FinancialMonteCarlo.jl simulates the underlying of a derivative under a stochastic model and
estimates prices, sensitivities and statistics of the payoff. Everything is built around a small
set of composable objects, so any model can be combined with any payoff, simulation method or
parallel mode, and parameters can be dual numbers to get derivatives by automatic differentiation.</p>
</div>
<nav class="pkg-links">
<a class="pkg-btn" href="https://github.com/rcalxrc08/FinancialMonteCarlo.jl"><svg viewBox="0 0 16 16" aria-hidden="true"><path d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013 8.013 0 0016 8c0-4.42-3.58-8-8-8z"/></svg> GitHub</a>
<a class="pkg-btn" href="https://github.com/rcalxrc08/FinancialMonteCarlo.jl/stargazers">★ Star</a>
<a class="pkg-btn" href="#Installation">Get started ↓</a>
<a class="pkg-btn" href="https://rcalxrc08.github.io/">All projects</a>
</nav>
```

```julia
using FinancialMonteCarlo

model = BlackScholesProcess(0.2, Underlying(100.0, 0.01))   # σ = 20%, S0 = 100, dividend 1%
rfCurve = ZeroRate(0.02)
mcConfig = MonteCarloConfiguration(10_000, 30)              # 10_000 paths, 30 time steps

pricer(model, rfCurve, mcConfig, EuropeanOption(1.0, 100.0))
```

## Features

- **Single name models**: Black-Scholes, Heston, Kou, Merton, Variance Gamma, Normal Inverse
  Gaussian, LogNormal mixture and shifted LogNormal mixture, plus the building block engines
  (Brownian motion, geometric Brownian motion, subordinated Brownian motion).
- **Multivariate models** through Gaussian copulas.
- **Time dependent zero rates and dividends**.
- **Payoffs**: forwards, European, binary, barrier (single and double), Asian (fixed and floating strike),
  American and Bermudan (Longstaff-Schwartz), basket European options and option smiles.
- **Contract algebra**: payoffs on the same underlying can be added and scaled into positions;
  portfolios on many underlyings are priced against a market data set.
- **Metrics**: price, variance, confidence interval, delta, rho, terminal distribution.
- **Simulation methods**: standard Monte Carlo, antithetic variates, control variates, Sobol
  sequences, any `AbstractRNG` (including [VectorizedRNG.jl](https://github.com/JuliaSIMD/VectorizedRNG.jl)).
- **Parallelism**: multi-threading, multi-processing ([Distributed](https://docs.julialang.org/en/v1/stdlib/Distributed/))
  and partial GPU support through [CUDA.jl](https://github.com/JuliaGPU/CUDA.jl).
- **Automatic differentiation** with [DualNumbers.jl](https://github.com/JuliaDiff/DualNumbers.jl),
  [HyperDualNumbers.jl](https://github.com/JuliaDiff/HyperDualNumbers.jl),
  [TaylorSeries.jl](https://github.com/JuliaDiff/TaylorSeries.jl),
  [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl) and
  [ReverseDiff.jl](https://github.com/JuliaDiff/ReverseDiff.jl).
- **Custom SDEs** defined with [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl).

## Installation

```julia
] add FinancialMonteCarlo
```

and, to run the test suite,

```julia
] test FinancialMonteCarlo
```

## Where to go next

- New to the package? Start with [Getting started](tutorials/getting_started.md).
- The [Manual](manual/overview.md) explains each building block in detail.
- [Extending the package](extending.md) shows how to add your own models and payoffs.
- The [API reference](api.md) lists every documented type and function.

!!! note "Keep in mind"
    First order automatic differentiation works with every payoff, including those for which the
    derivative of the Monte Carlo estimator is meaningless (e.g. binary and barrier options).
    Second order automatic differentiation works too, but the results are seldom useful.
    Multi-currency is not supported.
