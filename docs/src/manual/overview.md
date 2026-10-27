# Overview

Every computation in FinancialMonteCarlo.jl combines four objects:

| Object | Role | Examples |
|:--|:--|:--|
| model | how the underlying evolves | [`BlackScholesProcess`](@ref), [`HestonProcess`](@ref), [`GaussianCopulaNVariateProcess`](@ref) |
| zero rate curve | discounting and risk neutral drift | [`ZeroRate`](@ref), [`FinancialMonteCarlo.ImpliedZeroRate`](@ref) |
| Monte Carlo configuration | number of paths and steps, method, rng, parallelism | [`MonteCarloConfiguration`](@ref) |
| payoff | the contract | [`EuropeanOption`](@ref), [`BarrierOptionDownOut`](@ref), positions and portfolios |

and a **metric** applied to them, e.g. [`pricer`](@ref):

```julia
price = pricer(model, rfCurve, mcConfig, payoff)
```

## What happens inside `pricer`

1. the random number generator in `mcConfig` is reset to its seed, so results are reproducible;
2. `simulate` generates an `Nsim × (Nstep + 1)` matrix of paths of the underlying, from time
   zero up to the maturity of the payoff (or the largest maturity, when pricing several payoffs);
3. `payoff` maps each path to its discounted payout;
4. the average over the paths is returned.

All the other metrics ([`variance`](@ref), [`confinter`](@ref), [`delta`](@ref), [`rho`](@ref), ...)
reuse the same steps.

## Choosing the implementation by dispatch

There is no switch statement anywhere: the right simulation is selected by multiple dispatch on
the types of the arguments. For example, the same `simulate` call on a [`BlackScholesProcess`](@ref) runs

- standard Monte Carlo when `mcConfig` uses `StandardMC()`,
- antithetic sampling when it uses `AntitheticMC()`,
- Sobol low discrepancy sequences when it uses `SobolMode()`,
- a GPU kernel when its parallel mode is a `CudaMode` (with CUDA.jl loaded),

and accounts for time dependent rates when `rfCurve` is a curve instead of a scalar.

## Generic number types

All structs are parametric in the type of each of their numeric fields, and the numeric type of
the output is inferred from all the inputs. A dual number anywhere (volatility, spot, rate,
maturity, strike, ...) therefore propagates through the whole computation, see
[Automatic differentiation](@ref autodiff).

## Optional integrations

Some features are loaded only when the corresponding package is loaded in the same session
(through [Requires.jl](https://github.com/JuliaPackaging/Requires.jl)):

| Package | Enables |
|:--|:--|
| `CUDA` | `CudaMode` and GPU simulation of several models |
| `Distributed` | `MultiProcess` parallel mode |
| `DifferentialEquations` | `MonteCarloDiffEqModel`, see [DifferentialEquations.jl](@ref) |
| `DualNumbers`, `TaylorSeries` | compatibility methods for these number types |
| `VectorizedRNG` | fast vectorized random number generation |
