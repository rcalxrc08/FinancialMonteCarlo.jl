# FinancialMonteCarlo.jl

*Monte Carlo pricing and risk management of equity derivatives in Julia.*

```@raw html
<img src="assets/logo.png" alt="FinancialMonteCarlo.jl logo" width="120" style="float: right; margin-left: 1em;">
```

FinancialMonteCarlo.jl simulates the underlying of a derivative under a stochastic model and
estimates prices, sensitivities and statistics of the payoff. Everything is built around a small
set of composable objects, so any model can be combined with any payoff, simulation method or
parallel mode, and parameters can be dual numbers to get derivatives by automatic differentiation.

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
