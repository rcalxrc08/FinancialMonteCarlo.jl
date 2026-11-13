```@raw html
<div class="fc-page-head">
<p class="fc-kicker"><span>Manual 04</span> Processes and engines</p>
```

# Models


```@raw html
<div class="fc-meta"><span>Black–Scholes</span><span>Heston</span><span>Lévy</span><span>copulas</span></div>
<nav class="fc-toc">
<a href="#Processes-and-engines"><b>01</b>Processes and engines</a>
<a href="#Single-name-models"><b>02</b>Single name models</a>
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
    │   └── LevyProcess
    │       ├── ItoProcess             BlackScholesProcess, HestonProcess, LogNormalMixture, ShiftedLogNormalMixture
    │       ├── FiniteActivityProcess  KouProcess, MertonProcess
    │       └── InfiniteActivityProcess VarianceGammaProcess, NormalInverseGaussianProcess
    └── VectorialMonteCarloProcess
        └── NDimensionalMonteCarloProcess GaussianCopulaNVariateProcess, GaussianCopulaNVariateLogProcess
</pre>
```

## Single name models

All single name processes take their parameters followed by the [`Underlying`](@ref):

| Model | Constructor | Parameters |
|:--|:--|:--|
| Black-Scholes | `BlackScholesProcess(σ, underlying)` | volatility |
| Heston | `HestonProcess(σ, σ₀, λ, κ, ρ, θ, underlying)` | vol of vol, initial vol, market price of vol risk, mean reversion speed, correlation, long term variance |
| Kou | `KouProcess(σ, λ, p, λ₊, λ₋, underlying)` | volatility, jump intensity, probability of an up jump, rates of the exponential up and down jump sizes |
| Merton | `MertonProcess(σ, λ, μ_jump, σ_jump, underlying)` | volatility, jump intensity, mean and standard deviation of the normal jump sizes |
| Variance Gamma | `VarianceGammaProcess(σ, θ, κ, underlying)` | volatility, drift and variance rate of the gamma subordinator |
| Normal Inverse Gaussian | `NormalInverseGaussianProcess(σ, θ, κ, underlying)` | volatility, drift and variance rate of the inverse gaussian subordinator |
| LogNormal mixture | `LogNormalMixture(η, λ, underlying)` | vector of volatilities, vector of the first `length(η) - 1` weights |
| Shifted LogNormal mixture | `ShiftedLogNormalMixture(η, λ, α, underlying)` | as above, plus a shift |

Constructors check the parameters (positive volatilities and intensities, valid probabilities
and correlations, existence of the risk neutral drift, ...) and throw an `AssertionError` when they
are not admissible.

The Heston variance follows an Euler scheme floored at a tiny positive value (`sqrt(floatmin)` of
its number type, about `1e-154` in `Float64` and `1e-19` in `Float32`) instead of zero, where the
derivative of its square root is infinite. When the Feller condition `2κθ ≥ σ²` does not hold the
floor is reached often; it is far enough from the subnormal numbers that this costs nothing.

## Engines

| Engine | Constructor |
|:--|:--|
| Brownian motion with drift | `BrownianMotion(σ, μ)` |
| Geometric Brownian motion | `GeometricBrownianMotion(σ, μ, x0)` |
| Subordinated Brownian motion | `SubordinatedBrownianMotion(σ, μ, θ, subordinator::Distribution)` |

The drift `μ` can also be a curve, see [Market data](@ref).

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
