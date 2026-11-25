# FinancialMonteCarlo.jl <img src="etc/logo.png" width="40">

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://rcalxrc08.github.io/FinancialMonteCarlo.jl/stable/)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://rcalxrc08.github.io/FinancialMonteCarlo.jl/)
[![CI](https://github.com/rcalxrc08/FinancialMonteCarlo.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/rcalxrc08/FinancialMonteCarlo.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/rcalxrc08/FinancialMonteCarlo.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/rcalxrc08/FinancialMonteCarlo.jl)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

Monte Carlo pricing and risk of equity derivatives in Julia: single-name models from Black–Scholes to stochastic
volatility and jumps, Gaussian copulas for baskets, every payoff from forwards to American options, on threads,
processes or a GPU, and differentiable end to end by seven automatic differentiation backends.

Every metric is one function of four arguments, **model, zero rate curve, method, payoff**:

```julia
using FinancialMonteCarlo

model    = BlackScholesProcess(0.2)                             # σ: a model holds only its parameters
und      = Underlying(100.0, 0.01)                              # spot, dividend yield
rfCurve  = ZeroRate(0.02)
mcConfig = MonteCarloConfiguration(10_000, 30)                  # paths, steps
option   = EuropeanOption(1.0, 100.0)                           # maturity, strike

pricer(model, und, rfCurve, mcConfig, option)
confinter(model, und, rfCurve, mcConfig, option, 0.95)          # 95% confidence interval
```

## Installation

```julia
] add FinancialMonteCarlo
```

## What is inside

| | |
|:--|:--|
| **Models** | Black–Scholes, Bachelier, Heston, Kou, Merton, Variance Gamma, Normal Inverse Gaussian, lognormal mixture (Brigo–Mercurio), any of them shifted (`ShiftedProcess`); any SDE written with [StochasticDiffEq.jl](https://github.com/SciML/StochasticDiffEq.jl) |
| **Multivariate** | Gaussian copulas over any single-name models (`GaussianCopulaNVariateProcess`), basket payoffs |
| **Payoffs** | forwards, European and binary European, barriers, Asian, Bermudan and American (Longstaff–Schwartz), baskets; smiles of many strikes on one simulation |
| **Contracts** | positions (`2 * call - put`) and portfolios on many underlyings, priced against a market data set |
| **Metrics** | `pricer`, `variance`, `confinter`, `delta`, `rho`, `distribution`, `simulate` |
| **Sampling** | standard, antithetic, control variates, Sobol sequences, any `AbstractRNG` |
| **Execution** | serial, `MultiThreading`, `MultiProcess` and `CudaMode` (CUDA.jl): results do not depend on the number of threads or workers |
| **Market data** | constant or time dependent zero rates and dividends (`Curve`, `ImpliedCurve` from quoted zero rates) |
| **Number types** | the inputs decide: `Float32` stays `Float32`, dual numbers and Taylor series go through |

## Greeks by automatic differentiation

The price is a differentiable function of the inputs, so pathwise Greeks of any payoff come from the AD backend of
your choice:

```julia
using ForwardDiff
price(x) = pricer(BlackScholesProcess(x[2]), Underlying(x[1], 0.01), ZeroRate(x[3]), mcConfig, option)
ForwardDiff.gradient(price, [100.0, 0.2, 0.02])   # delta, vega, rho
```

| Backend | Simulations | Notes |
|:--|:--|:--|
| [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl) | ✓ CPU and GPU | |
| [Enzyme.jl](https://github.com/EnzymeAD/Enzyme.jl) | ✓ CPU, forward and reverse | closed-form adjoints of the inner kernels; no runtime activity needed, curves included |
| [Mooncake.jl](https://github.com/chalk-lab/Mooncake.jl) | ✓ CPU, forward and reverse | closed-form adjoints of the inner kernels |
| [ReverseDiff.jl](https://github.com/JuliaDiff/ReverseDiff.jl) | ✓ CPU | |
| [DualNumbers.jl](https://github.com/JuliaDiff/DualNumbers.jl), [HyperDualNumbers.jl](https://github.com/JuliaDiff/HyperDualNumbers.jl), [TaylorSeries.jl](https://github.com/JuliaDiff/TaylorSeries.jl) | ✓ CPU | second order through hyper-duals and Taylor series |
| [Zygote.jl](https://github.com/FluxML/Zygote.jl) | not the simulations | the models, `set_parameters!`, the curves and their integrals, with reverse rules: what a calibration loss differentiates |

First order AD runs for every payoff, also where the pathwise derivative is meaningless (binary and barrier options);
second order results are mostly useless. See the [autodiff manual](https://rcalxrc08.github.io/FinancialMonteCarlo.jl/stable/manual/autodiff/)
for timings and caveats.

## Positions, smiles and portfolios

```julia
position = EuropeanOption(1.0, 100.0) - 2 * BinaryEuropeanOption(1.0, 110.0)
pricer(model, und, rfCurve, mcConfig, position)

smile = FinancialMonteCarlo.EuropeanOptionSmile(1.0, [90.0, 100.0, 110.0], [false, true, true])
pricer(model, und, rfCurve, mcConfig, smile)       # one price per strike, from one simulation

market    = ("ENI" → (model, und)) + ("AAPL" → (HestonProcess(0.3, 0.2, 0.0, 1.5, -0.5, 0.04), Underlying(150.0)))
portfolio = ("ENI" → position) + ("AAPL" → 3.0 * AsianFixedStrikeOption(1.0, 150.0))
pricer(market, rfCurve, mcConfig, portfolio)       # → is typed \to<TAB>
```

## Curves

Rates and dividends can depend on time. A curve holds instantaneous rates, linear between its nodes and flat outside;
`ImpliedCurve` builds it from zero rates quoted at given maturities. `ZeroRate` and `Underlying` take a number or a
curve:

```julia
rfCurve   = ZeroRate(FinancialMonteCarlo.ImpliedCurve([0.01, 0.015, 0.02], [0.5, 1.0, 2.0]))
dividends = FinancialMonteCarlo.ImpliedCurve([0.010, 0.012, 0.011], [0.5, 1.0, 2.0])
pricer(model, Underlying(100.0, dividends), rfCurve, mcConfig, option)
```

The nodes are static vectors, so an underlying on a curve is an `isbits` value: Enzyme differentiates it without runtime
activity, and a GPU kernel can capture it. There is no multi-currency support.

## On the GPU

```julia
using CUDA
gpuConfig = MonteCarloConfiguration(100_000, 50; mode = FinancialMonteCarlo.CudaMode())
pricer(BlackScholesProcess(0.2f0), Underlying(100.0f0, 0.01f0), ZeroRate(0.02f0), gpuConfig, EuropeanOption(1.0f0, 100.0f0))
```

The paths are simulated and the payoffs evaluated on the device; in `Float32` a laptop GPU runs 100 000 paths × 50
steps in about 6 ms. See the [performance manual](https://rcalxrc08.github.io/FinancialMonteCarlo.jl/stable/manual/performance/).

## The ecosystem

The `pricer` interface is shared by the other packages of the author, which add their own methods on the same models
and payoffs:

- [FinancialFFT.jl](https://github.com/rcalxrc08/FinancialFFT.jl): Fourier methods (Lewis, COS, Carr–Madan) and
  densities of the log-return;
- [VibratoMonteCarlo.jl](https://github.com/rcalxrc08/VibratoMonteCarlo.jl): smooth Greeks of discontinuous payoffs
  (Vibrato, likelihood ratio method);
- [FinancialCalibrator.jl](https://github.com/rcalxrc08/FinancialCalibrator.jl): calibration to implied volatility
  surfaces, differentiable in the quotes;
- [FinancialToolbox.jl](https://github.com/rcalxrc08/FinancialToolbox.jl): Black–Scholes closed forms and implied
  volatility, the references of the test suite.

## Tests

```julia
] test FinancialMonteCarlo
```
