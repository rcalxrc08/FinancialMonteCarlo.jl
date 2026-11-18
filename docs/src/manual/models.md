```@raw html
<div class="fc-page-head">
<p class="fc-kicker"><span>Manual 04</span> How models are built and simulated</p>
```

# Models

How a model is put together and simulated: processes and engines, the constructors, `simulate`, multivariate
models and parameters. Each model, with its dynamics and what it supports, has its card in the
[Model bestiary](@ref).

```@raw html
<div class="fc-meta"><span>processes</span><span>engines</span><span>simulate</span><span>parameters</span></div>
<nav class="fc-toc">
<a href="#Processes-and-engines"><b>01</b>Processes and engines</a>
<a href="#Constructors"><b>02</b>Constructors</a>
<a href="#Engines"><b>03</b>Engines</a>
<a href="#Simulation"><b>04</b>Simulation</a>
<a href="#Multivariate-models"><b>05</b>Multivariate models</a>
<a href="#Parameters"><b>06</b>Parameters</a>
</nav>
</div>
```
## Processes and engines

Models come in two flavours:

- **processes** (subtypes of `AbstractMonteCarloProcess`) describe a financial asset: they hold an
  [`Underlying`](@ref) and, combined with a zero rate curve, simulate the asset under the risk
  neutral measure. They are what you pass to the metrics;
- **engines** (subtypes of `AbstractMonteCarloEngine`) are plain stochastic processes, with no
  notion of rates or dividends. Processes are built on top of them: e.g. a
  [`KouProcess`](@ref) computes its risk neutral drift, simulates a [`BrownianMotion`](@ref),
  adds the jumps and exponentiates.

```@raw html
<pre>
BaseProcess
├── AbstractMonteCarloEngine           BrownianMotion, GeometricBrownianMotion, SubordinatedBrownianMotion
└── AbstractMonteCarloProcess
    ├── ScalarMonteCarloProcess
    │   ├── ShiftedProcess
    │   └── LevyProcess
    │       ├── ItoProcess             BlackScholesProcess, BachelierProcess, HestonProcess, LogNormalMixture
    │       ├── FiniteActivityProcess  KouProcess, MertonProcess
    │       └── InfiniteActivityProcess VarianceGammaProcess, NormalInverseGaussianProcess
    └── VectorialMonteCarloProcess
        └── NDimensionalMonteCarloProcess GaussianCopulaNVariateProcess, GaussianCopulaNVariateLogProcess
</pre>
```

## Constructors

All single name processes take their parameters followed by the [`Underlying`](@ref), e.g.
`HestonProcess(σ, σ₀, λ, κ, ρ, θ, underlying)`; the [Model bestiary](@ref) lists them all. A wrapper such as
[`ShiftedProcess`](@ref) takes a model instead, and the copulas take several.

Constructors check the parameters (positive volatilities and intensities, valid probabilities
and correlations, existence of the risk neutral drift, ...) and throw a `DomainError` when they
are not admissible (a `DimensionMismatch` for a wrong number of parameters or mismatched sizes, an `ArgumentError`
for the other invalid inputs), so a caller can catch exactly these and let a bug propagate. They are fully parametric in the number type of each parameter, so a model holds dual numbers,
`Float32` or any other number type, and the paths and prices follow it.

## Engines

The engines are plain stochastic processes, with no rates or dividends: [`simulate`](@ref) takes no
zero rate curve for them. The models are built on top of them, and so can new ones.

| Engine | Dynamics | Used by |
|:--|:--|:--|
| [`BrownianMotion`](@ref)`(σ, μ)` | ``X_t = μ t + σ W_t`` | Bachelier; Merton, Kou (plus their jumps) |
| [`GeometricBrownianMotion`](@ref)`(σ, μ, x0)` | ``dX_t = μ X_t\,dt + σ X_t\,dW_t`` | Black-Scholes |
| [`SubordinatedBrownianMotion`](@ref)`(σ, drift, θ, subordinator)` | ``X_{t+Δt} = X_t + \text{drift}\,Δt + θ\,Δτ + σ\sqrt{Δτ}\,Z`` | Variance Gamma, Normal Inverse Gaussian |

The drift can also be a curve, see [Market data](@ref). The `subordinator` of a [`SubordinatedBrownianMotion`](@ref)
is any continuous univariate `Distribution` of the increment over one step, so other time changed Lévy models need
only their law and their martingale drift. Heston and the lognormal mixture have Euler schemes of their own, on the
hooks of the package (antithetic sides, draws on the device of the mode), see [Extending the package](@ref).

## Simulation

[`simulate`](@ref) returns the simulated paths, one per row:

```julia
S = simulate(model, rfCurve, mcConfig, T)        # processes: Nsim × (Nstep + 1) matrix
X = simulate(BrownianMotion(0.2, 0.0), mcConfig, T)  # engines need no zero rate curve
```

Multivariate processes return a vector of such matrices, one per component.

## Multivariate models

[`GaussianCopulaNVariateProcess`](@ref) and [`GaussianCopulaNVariateLogProcess`](@ref) couple
single name models through a Gaussian copula with correlation matrix `ρ`:

```julia
GaussianCopulaNVariateProcess(ρ, model_1, model_2, model_3)
GaussianCopulaNVariateProcess(model_1, model_2, ρ₁₂)   # bivariate shortcut
```

The log variant applies the copula to the log returns of the components. See
[Portfolios and multivariate](@ref) for how to price contracts on them.

## Parameters

[`FinancialMonteCarlo.get_parameters`](@ref) returns the parameters of a model (or of a payoff) as
a vector, and [`FinancialMonteCarlo.set_parameters!`](@ref) builds a new model with the given
parameters. Together they make calibration loops independent of the model. Both are type stable:
their result is inferred from the type of the model, also when the parameters are dual numbers, so
the code that follows them in a loss function is inferred too. The parameters are the fields of the
model in their order, without the underlying (and without the models of a copula).

See [Models](@ref api-models) in the API reference.
