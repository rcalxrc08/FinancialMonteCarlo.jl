```@raw html
<div class="fc-page-head">
<p class="fc-kicker"><span>Manual 05</span> Every model, one card each</p>
```

# Model bestiary

Every model the package simulates: its dynamics under the risk neutral measure, its constructor and
parameters, how it is simulated and what it supports. [Models](@ref) explains how processes and
engines fit together; this page is the catalogue.

```@raw html
<div class="fc-meta"><span>diffusions</span><span>stochastic volatility</span><span>jumps</span><span>Lévy</span><span>mixtures</span><span>copulas</span></div>
<nav class="fc-toc">
<a href="#At-a-glance"><b>01</b>At a glance</a>
<a href="#Black-Scholes"><b>02</b>Black-Scholes</a>
<a href="#Heston"><b>03</b>Heston</a>
<a href="#Merton"><b>04</b>Merton</a>
<a href="#Kou"><b>05</b>Kou</a>
<a href="#Variance-Gamma"><b>06</b>Variance Gamma</a>
<a href="#Normal-Inverse-Gaussian"><b>07</b>Normal Inverse Gaussian</a>
<a href="#LogNormal-mixture"><b>08</b>LogNormal mixture</a>
<a href="#Shifted-LogNormal-mixture"><b>09</b>Shifted LogNormal mixture</a>
<a href="#Gaussian-copulas"><b>10</b>Gaussian copulas</a>
<a href="#Engines"><b>11</b>Engines</a>
<a href="#Any-SDE"><b>12</b>Any SDE</a>
<a href="#Not-yet-here"><b>13</b>Not yet here</a>
</nav>
</div>
```

## At a glance

In every model ``r_t`` is the zero rate and ``d_t`` the dividend yield of the [`Underlying`](@ref),
both constants or curves (see [Market data](@ref)), and the drift makes the discounted price, with
dividends, a martingale. Every single name model works with the three samplings (standard,
antithetic, Sobol), with rate and dividend curves and in the batched modes (`MultiThreading`,
`MultiProcess`).

| Model | Family | Parameters | GPU | Enzyme | Closed form in FinancialFFT |
|:--|:--|:--|:--:|:--:|:--:|
| [`BlackScholesProcess`](@ref) | diffusion | 1 | tested | yes | yes |
| [`HestonProcess`](@ref) | stochastic volatility | 6 | | yes | yes |
| [`MertonProcess`](@ref) | jump diffusion | 4 | tested | yes | yes |
| [`KouProcess`](@ref) | jump diffusion | 5 | tested | yes | yes |
| [`VarianceGammaProcess`](@ref) | infinite activity Lévy | 3 | tested | no | yes |
| [`NormalInverseGaussianProcess`](@ref) | infinite activity Lévy | 3 | tested | yes | yes |
| [`LogNormalMixture`](@ref) | mixture | ``2n - 1`` | | yes | |
| [`ShiftedLogNormalMixture`](@ref) | mixture | ``2n`` | | | |
| [`GaussianCopulaNVariateProcess`](@ref) | multivariate | ``ρ`` + marginals | | | |
| [`GaussianCopulaNVariateLogProcess`](@ref) | multivariate | ``ρ`` + marginals | | | |

- **GPU**: covered by the suite in `test/cuda/` (run by hand, it needs a GPU); the other models may
  work in `CudaMode` but are not tested there. See [Performance](@ref).
- **Enzyme**: forward and reverse mode agree with ForwardDiff (see [Automatic differentiation](@ref autodiff)).
  ForwardDiff, ReverseDiff, Mooncake, DualNumbers and TaylorSeries work with every model.
- **Closed form**: [FinancialFFT.jl](https://github.com/rcalxrc08/FinancialFFT.jl) prices European
  options from the characteristic function of the model, which gives a reference price for the
  Monte Carlo one (and for Black-Scholes, `blsprice` of
  [FinancialToolbox.jl](https://github.com/rcalxrc08/FinancialToolbox.jl)).

## Black-Scholes

```math
dS_t = (r_t - d_t) S_t\,dt + σ S_t\,dW_t
```

```julia
BlackScholesProcess(σ, underlying)
BlackScholesProcess(0.2, Underlying(100.0, 0.01))
```

| Parameter | Meaning | Constraint |
|:--|:--|:--|
| `σ` | volatility | ``σ > 0`` |

Simulated exactly (a [`GeometricBrownianMotion`](@ref) with the risk neutral drift), so a European
payoff needs one time step. The reference for everything else: its European, binary and implied
volatility formulas are in FinancialToolbox (`blsprice`, `blsbin`, `blsimpv`).

## Heston

```math
dS_t = (r_t - d_t) S_t\,dt + \sqrt{v_t}\,S_t\,dW_t, \qquad
dv_t = κ^*(θ^* - v_t)\,dt + σ \sqrt{v_t}\,dZ_t, \qquad
d\langle W, Z\rangle_t = ρ\,dt
```

with ``v_0 = σ_0^2``, ``κ^* = κ + λ`` and ``θ^* = κθ / (κ + λ)``.

```julia
HestonProcess(σ, σ₀, λ, κ, ρ, θ, underlying)
HestonProcess(0.5, 0.2, 0.0, 2.0, -0.7, 0.04, Underlying(100.0, 0.01))
```

| Parameter | Meaning | Constraint |
|:--|:--|:--|
| `σ` | volatility of the variance | ``σ > 0`` |
| `σ₀` | initial volatility | ``σ_0 > 0`` |
| `λ` | market price of volatility risk | ``κ + λ ≠ 0`` |
| `κ` | speed of mean reversion | |
| `ρ` | correlation of the two Brownian motions | ``-1 ≤ ρ ≤ 1`` |
| `θ` | long term variance | |

Euler scheme on the log price and on the variance, which is floored at a tiny positive value
(``\sqrt{\text{floatmin}}`` of its number type) instead of zero. A smile with ``ρ < 0`` slopes down.
When the Feller condition ``2κθ ≥ σ^2`` is badly violated the Euler scheme is biased: see
[#55](https://github.com/rcalxrc08/FinancialMonteCarlo.jl/issues/55), and use more time steps.

## Merton

```math
\log \frac{S_t}{S_0} = \int_0^t (r_s - d_s)\,ds + ω t + σ W_t + \sum_{i=1}^{N_t} Y_i, \qquad
Y_i \sim \mathcal{N}(μ_J, σ_J^2)
```

with ``N_t`` a Poisson process of intensity ``λ`` and ``ω = -σ^2/2 - λ\left(e^{μ_J + σ_J^2/2} - 1\right)``.

```julia
MertonProcess(σ, λ, μ_jump, σ_jump, underlying)
MertonProcess(0.2, 3.0, -0.05, 0.1, Underlying(100.0, 0.01))
```

| Parameter | Meaning | Constraint |
|:--|:--|:--|
| `σ` | volatility of the diffusion | ``σ > 0`` |
| `λ` | jump intensity (jumps per year) | ``λ > 0`` |
| `μ_jump` | mean of the jumps of the log price | |
| `σ_jump` | standard deviation of the jumps of the log price | ``σ_J > 0`` |

A finite activity Lévy process: a Brownian motion plus a compound Poisson process
(`FiniteActivityProcess`, [fa_levy.jl](https://github.com/rcalxrc08/FinancialMonteCarlo.jl/blob/master/src/models/fa_levy.jl)).
The jumps are simulated exactly, not discretised; with antithetic and Sobol sampling their count,
times and sizes are drawn from uniforms by inversion, so that the paths of a pair mirror each other
and the Sobol points cover the jumps too.

## Kou

```math
\log \frac{S_t}{S_0} = \int_0^t (r_s - d_s)\,ds + ω t + σ W_t + \sum_{i=1}^{N_t} Y_i, \qquad
Y_i = \begin{cases} +\text{Exp}(λ_+) & \text{with probability } p \\ -\text{Exp}(λ_-) & \text{with probability } 1 - p \end{cases}
```

with ``N_t`` a Poisson process of intensity ``λ`` and
``ω = -σ^2/2 - λ\left(\frac{p λ_+}{λ_+ - 1} + \frac{(1 - p) λ_-}{λ_- + 1} - 1\right)``.

```julia
KouProcess(σ, λ, p, λ₊, λ₋, underlying)
KouProcess(0.2, 5.0, 0.3, 30.0, 20.0, Underlying(100.0, 0.01))
```

| Parameter | Meaning | Constraint |
|:--|:--|:--|
| `σ` | volatility of the diffusion | ``σ > 0`` |
| `λ` | jump intensity | ``λ > 0`` |
| `p` | probability of an up jump | ``0 ≤ p ≤ 1`` |
| `λ₊` | rate of the exponential up jumps | ``λ_+ > 1`` (finite forward) |
| `λ₋` | rate of the exponential down jumps | ``λ_- > 0`` |

Simulated like Merton, by the same finite activity code: only the law of the jump sizes differs.

## Variance Gamma

```math
\log \frac{S_t}{S_0} = \int_0^t (r_s - d_s)\,ds + ω t + θ G_t + σ W_{G_t}, \qquad
G_t \sim \text{Gamma}\left(\frac{t}{κ}, κ\right)
```

with ``ω = \frac{1}{κ} \log\left(1 - \frac{σ^2 κ}{2} - θκ\right)``: a Brownian motion with drift
evaluated at a gamma subordinator of unit mean rate and variance rate ``κ``.

```julia
VarianceGammaProcess(σ, θ, κ, underlying)
VarianceGammaProcess(0.2, -0.1, 0.2, Underlying(100.0, 0.01))
```

| Parameter | Meaning | Constraint |
|:--|:--|:--|
| `σ` | volatility of the subordinated Brownian motion | ``σ > 0`` |
| `θ` | drift of the subordinated Brownian motion (skew) | ``1 - σ^2 κ / 2 - θ κ > 0`` |
| `κ` | variance rate of the subordinator (kurtosis) | ``κ > 0`` |

An infinite activity Lévy process (`InfiniteActivityProcess`,
[ia_levy.jl](https://github.com/rcalxrc08/FinancialMonteCarlo.jl/blob/master/src/models/ia_levy.jl)),
exact at the time steps: a [`SubordinatedBrownianMotion`](@ref) whose subordinator increments are
gamma draws (on the GPU, by a Marsaglia-Tsang kernel). Enzyme cannot differentiate it, because the
gamma sampler of Distributions.jl draws scalar exponentials.

## Normal Inverse Gaussian

```math
\log \frac{S_t}{S_0} = \int_0^t (r_s - d_s)\,ds + ω t + θ I_t + σ W_{I_t}, \qquad
I_t \sim \text{InverseGaussian}\left(t, \frac{t^2}{κ}\right)
```

with ``ω = \frac{1}{κ}\left(\sqrt{1 - σ^2 κ - 2θκ} - 1\right)``: the same construction as Variance
Gamma with an inverse Gaussian subordinator.

```julia
NormalInverseGaussianProcess(σ, θ, κ, underlying)
NormalInverseGaussianProcess(0.2, -0.1, 0.2, Underlying(100.0, 0.01))
```

| Parameter | Meaning | Constraint |
|:--|:--|:--|
| `σ` | volatility of the subordinated Brownian motion | ``σ > 0`` |
| `θ` | drift of the subordinated Brownian motion | ``1 - σ^2 κ - 2θκ > 0`` |
| `κ` | variance rate of the subordinator | ``κ > 0`` |

On the GPU the subordinator is drawn by a Michael-Schucany-Haas kernel.

## LogNormal mixture

```math
S_t = \sum_{i=1}^n λ_i S^i_t, \qquad dS^i_t = (r_t - d_t) S^i_t\,dt + η_i S^i_t\,dW^i_t, \qquad S^i_0 = S_0
```

```julia
LogNormalMixture(η, λ, underlying)
LogNormalMixture([0.2, 0.3], [0.4], Underlying(100.0, 0.01))   # weights 0.4 and 0.6
```

| Parameter | Meaning | Constraint |
|:--|:--|:--|
| `η` | vector of the ``n`` volatilities | all positive |
| `λ` | vector of the first ``n - 1`` weights; the last is ``1 - \sum λ`` | positive, ``\sum λ ≤ 1`` |

!!! warning
    The paths are a weighted sum of geometric Brownian motions, which is not the lognormal mixture
    of Brigo and Mercurio (whose density is the mixture of the lognormal densities, and whose
    European price is the weighted sum of Black-Scholes prices): see
    [#56](https://github.com/rcalxrc08/FinancialMonteCarlo.jl/issues/56).

## Shifted LogNormal mixture

```math
S_t = \sum_{i=1}^n λ_i S^i_t + α S_0\, e^{\int_0^t (r_s - d_s)\,ds}, \qquad S^i_0 = (1 - α) S_0
```

with the ``S^i`` of the [LogNormal mixture](@ref).

```julia
ShiftedLogNormalMixture(η, λ, α, underlying)
ShiftedLogNormalMixture([0.2, 0.3], [0.4], 0.1, Underlying(100.0, 0.01))
```

| Parameter | Meaning | Constraint |
|:--|:--|:--|
| `η`, `λ` | as in the [LogNormal mixture](@ref) | |
| `α` | shift, as a fraction of the spot | |

The deterministic shift grows at the risk neutral rate, so the forward is the one of the other
models. It shares the caveat of the LogNormal mixture.

## Gaussian copulas

```julia
GaussianCopulaNVariateProcess(ρ, model_1, model_2, model_3)
GaussianCopulaNVariateProcess(model_1, model_2, ρ₁₂)     # bivariate shortcut
GaussianCopulaNVariateLogProcess(ρ, model_1, model_2)    # copula on the log returns
```

| Parameter | Meaning | Constraint |
|:--|:--|:--|
| `ρ` | correlation matrix, one row per component | symmetric positive definite |
| `models` | any single name models above, also different ones | |

Each component is simulated with its own model, so it keeps its marginal law; at every time step
its values are then replaced by their empirical quantiles at the uniforms of a Gaussian copula with
correlation ``ρ``. They price baskets
([`EuropeanOptionND`](@ref)) and the positions of a [Portfolios and multivariate](@ref) set.
Known limitations: the ranks are drawn again at each step, so a path does not stay continuous
([#23](https://github.com/rcalxrc08/FinancialMonteCarlo.jl/issues/23)), and the reordering sorts
every column at every step ([#46](https://github.com/rcalxrc08/FinancialMonteCarlo.jl/issues/46)).

## Engines

The engines are plain stochastic processes, with no rates or dividends: [`simulate`](@ref) takes no
zero rate curve for them. Models are built on top of them, and so can new ones.

| Engine | Dynamics | Used by |
|:--|:--|:--|
| [`BrownianMotion`](@ref)`(σ, μ)` | ``X_t = μ t + σ W_t`` | Merton, Kou (plus their jumps) |
| [`GeometricBrownianMotion`](@ref)`(σ, μ, x0)` | ``dX_t = μ X_t\,dt + σ X_t\,dW_t`` | Black-Scholes, LogNormal mixtures |
| [`SubordinatedBrownianMotion`](@ref)`(σ, drift, θ, subordinator)` | ``X_{t+Δt} = X_t + \text{drift}\,Δt + θ\,Δτ + σ\sqrt{Δτ}\,Z`` | Variance Gamma, Normal Inverse Gaussian |

A drift can be a curve. The `subordinator` of a [`SubordinatedBrownianMotion`](@ref) is any
continuous univariate `Distribution` of the increment over one step, so other time changed Lévy
models need only their law and their martingale drift.

## Any SDE

With DifferentialEquations.jl loaded, `FinancialMonteCarlo.MonteCarloDiffEqModel` turns an
`EnsembleProblem` (an SDE, or a jump diffusion) into a model priced like the others: see
[DifferentialEquations.jl](@ref). It is the quickest way to try a model that is not here, at the
cost of the speed, of the GPU and of the samplings of the native models.

## Not yet here

Bachelier, SABR, CEV, Bates, displaced diffusions and other models that the current building blocks
make easy to add are collected in
[#62](https://github.com/rcalxrc08/FinancialMonteCarlo.jl/issues/62). To write one yourself, see
[Extending the package](@ref).
