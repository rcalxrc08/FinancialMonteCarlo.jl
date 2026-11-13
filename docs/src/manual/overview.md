```@raw html
<div class="fc-page-head">
<p class="fc-kicker"><span>Manual 01</span> A few composable objects</p>
```

# Overview

Every computation in FinancialMonteCarlo.jl combines four objects:

```@raw html
<div class="fc-meta"><span>model</span><span>curve</span><span>configuration</span><span>payoff</span></div>
<nav class="fc-toc">
<a href="#What-happens-inside-pricer"><b>01</b>What happens inside pricer</a>
<a href="#Choosing-the-implementation-by-dispatch"><b>02</b>Choosing the implementation by dispatch</a>
<a href="#Generic-number-types"><b>03</b>Generic number types</a>
<a href="#Optional-integrations"><b>04</b>Optional integrations</a>
</nav>
</div>
```

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

- standard Monte Carlo when `mcConfig` uses `StandardMC()` sampling,
- antithetic sampling when it uses `AntitheticMC()`,
- Sobol low discrepancy sequences when it uses `SobolMC()`,
- a GPU kernel when its mode is a `CudaMode` (with CUDA.jl loaded),

and accounts for time dependent rates when `rfCurve` is a curve instead of a scalar.

## Generic number types

All structs are parametric in the type of each of their numeric fields, and the numeric type of
the output is inferred from all the inputs. A dual number anywhere (volatility, spot, rate,
maturity, strike, ...) therefore propagates through the whole computation, see
[Automatic differentiation](@ref autodiff).

## Optional integrations

Some features are package extensions, loaded only when the corresponding package is loaded in the
same session:

| Package | Enables |
|:--|:--|
| `CUDA` | GPU simulation of several models with `CudaMode` |
| `Distributed` | `MultiProcess` parallel mode |
| `DifferentialEquations` | `MonteCarloDiffEqModel`, see [DifferentialEquations.jl](@ref) |
| `DualNumbers`, `TaylorSeries` | compatibility methods for these number types |
| `VectorizedRNG` | fast vectorized random number generation |
