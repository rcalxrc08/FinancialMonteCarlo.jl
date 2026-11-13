```@raw html
<div class="fc-hero">
<div class="fc-hero-text">
<p class="fc-eyebrow"><img src="assets/logo.png" alt="" width="28" height="28"> FinancialMonteCarlo.jl</p>
```

# Monte Carlo Pricing and Risk for Equity Derivatives

```@raw html
<p class="fc-lede">Simulate Black–Scholes, Heston and Lévy models, price anything from forwards to American options, and get sensitivities by automatic differentiation, on threads, processes or a GPU.</p>
<nav class="fc-cta">
<a class="fc-btn fc-btn-primary" href="#Quick-start">Get Started <span aria-hidden="true">→</span></a>
<a class="fc-btn" href="https://github.com/rcalxrc08/FinancialMonteCarlo.jl"><svg viewBox="0 0 16 16" aria-hidden="true"><path d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013 8.013 0 0016 8c0-4.42-3.58-8-8-8z"/></svg> View on GitHub</a>
</nav>
<pre class="fc-hero-code"><code><span class="k">model</span> = HestonProcess(0.4, 0.2, 0.0, 1.5, -0.7, 0.04, Underlying(100.0, 0.01))
<span class="k">mc</span>    = MonteCarloConfiguration(100_000, 50)
<span class="k">price</span> = pricer(model, ZeroRate(0.02), mc, AmericanOption(1.0, 100.0, false))</code></pre>
</div>
<figure class="fc-hero-figure">
<img src="assets/figures/hero-paths-dark.svg" alt="Simulated paths of a Heston model, coloured by their final value">
<img src="assets/figures/hero-paths-light.svg" alt="Simulated paths of a Heston model, coloured by their final value">
<figcaption><span class="fc-legend"><i class="fc-ramp"></i>coloured by the spot at one year</span>
160 paths of a Heston model over one year, simulated with FinancialMonteCarlo.jl when this site was built.</figcaption>
</figure>
</div>
```


```@raw html
<div class="fc-stats"><div class="fc-stat is-teal"><div class="fc-stat-value">6.2<span> ms</span></div><div class="fc-stat-label">100 000 paths × 50 steps on a laptop GPU, Float32</div></div><div class="fc-stat is-amber"><div class="fc-stat-value">39<span>×</span></div><div class="fc-stat-label">GPU in Float32 against one CPU thread, at a million paths</div></div><div class="fc-stat "><div class="fc-stat-value">9<span> models</span></div><div class="fc-stat-label">Black–Scholes, Heston, Kou, Merton, VG, NIG, mixtures, plus copulas</div></div><div class="fc-stat "><div class="fc-stat-value">7<span></span></div><div class="fc-stat-label">AD backends through the simulation, Enzyme and Mooncake included</div></div></div>
<p class="fc-machine">European call under Black–Scholes, Intel Core i7-8550U and NVIDIA MX150, see <a href="manual/performance/">Performance</a>.</p>
```

## What it does

```@raw html
<p class="fc-section-lede">Everything is built around a small set of composable objects: <em>any model with any payoff, simulation method or parallel mode</em>, and parameters that can be dual numbers.</p>
<div class="fc-grid fc-features">
<article class="fc-card"><div class="fc-icon"><svg viewBox="0 0 32 32"><path d="M4 16 9 13l4 5 5-9 5 4 5-6"/><path d="M4 16 9 19l4-2 5 6 5-3 5 5"/><circle cx="4" cy="16" r="2"/></svg></div>
<h3>Single-name models</h3><p>Black–Scholes, Heston, Kou, Merton, Variance Gamma, NIG, lognormal and shifted lognormal mixtures, on Brownian engines.</p><p class="fc-why">From textbook to jumps and stochastic volatility with one interface.</p></article>
<article class="fc-card"><div class="fc-icon"><svg viewBox="0 0 32 32"><circle cx="16" cy="6" r="2.5"/><circle cx="8" cy="18" r="2.5"/><circle cx="24" cy="18" r="2.5"/><path d="M16 8.5 8 15.5M16 8.5l8 7M8 20.5v6M24 20.5v6"/></svg></div>
<h3>Multivariate</h3><p>Gaussian copulas over any single-name models, for baskets and portfolios.</p><p class="fc-why">Correlated underlyings without a new model.</p></article>
<article class="fc-card"><div class="fc-icon"><svg viewBox="0 0 32 32"><rect x="4" y="4" width="7" height="7" rx="1.5"/><rect x="12.5" y="4" width="7" height="7" rx="1.5"/><rect x="21" y="4" width="7" height="7" rx="1.5"/><rect x="4" y="12.5" width="7" height="7" rx="1.5"/><rect x="12.5" y="12.5" width="7" height="7" rx="1.5" class="hl"/><rect x="21" y="12.5" width="7" height="7" rx="1.5"/><path d="M16 21v7"/></svg></div>
<h3>Payoffs</h3><p>Forwards, European, binary, barrier, Asian, Bermudan and American (Longstaff–Schwartz), baskets and smiles.</p><p class="fc-why">Every metric works for every payoff.</p></article>
<article class="fc-card"><div class="fc-icon"><svg viewBox="0 0 32 32"><path d="M4 9h24M4 16h24M4 23h24"/><circle cx="10" cy="9" r="2.4"/><circle cx="20" cy="16" r="2.4"/><circle cx="14" cy="23" r="2.4"/></svg></div>
<h3>Positions and portfolios</h3><p>Payoffs combine with <code>+</code> and <code>*</code>; portfolios on many underlyings price against a market data set.</p><p class="fc-why">A book priced in one call, from one stream of random numbers.</p></article>
<article class="fc-card"><div class="fc-icon"><svg viewBox="0 0 32 32"><path d="M6 26V14M12 26V8M18 26V17M24 26V11"/><path d="M4 27h24"/></svg></div>
<h3>Sampling</h3><p>Standard, antithetic, control variates, Sobol sequences, any <code>AbstractRNG</code>.</p><p class="fc-why">Variance reduction without touching the model.</p></article>
<article class="fc-card"><div class="fc-icon"><svg viewBox="0 0 32 32"><rect x="7" y="7" width="18" height="18" rx="3"/><rect x="12" y="12" width="8" height="8" rx="1"/><path d="M11 3v4M16 3v4M21 3v4M11 25v4M16 25v4M21 25v4M3 11h4M3 16h4M3 21h4M25 11h4M25 16h4M25 21h4"/></svg></div>
<h3>Threads, processes, GPU</h3><p><code>MultiThreading</code>, <code>MultiProcess</code> and <code>CudaMode</code>: the mode is part of the configuration.</p><p class="fc-why">Results that do not depend on the number of threads or workers.</p></article>
<article class="fc-card"><div class="fc-icon"><svg viewBox="0 0 32 32"><path d="M4 26 16 6l12 20z"/><path d="M10 26 16 16l6 10"/></svg></div>
<h3>Automatic differentiation</h3><p>DualNumbers, HyperDualNumbers, TaylorSeries, ForwardDiff, ReverseDiff, Mooncake, Enzyme.</p><p class="fc-why">Pathwise Greeks of any payoff from the price function.</p></article>
<article class="fc-card"><div class="fc-icon"><svg viewBox="0 0 32 32"><path d="M3 16c3-9 6-9 9 0s6 9 9 0 6-9 8-4"/><circle cx="12" cy="16" r="2" /></svg></div>
<h3>Custom SDEs</h3><p>Any model written with DifferentialEquations.jl, priced like the built-in ones.</p><p class="fc-why">Prototype a model before writing its simulation.</p></article>
</div>
```

## How it fits together

```@raw html
<p class="fc-section-lede">Every metric reduces to one pipeline: simulate the paths, evaluate the payoff on them, average.</p>
<div class="fc-arch">
<ol class="fc-flow" style="--fc-steps: 6">
<li><span class="fc-node-k">01</span><b>Market data</b><span><code>Underlying</code>, <code>ZeroRate</code>, curves</span></li>
<li><span class="fc-node-k">02</span><b>Model</b><span><code>HestonProcess</code>, <code>KouProcess</code>, …</span></li>
<li class="fc-core"><span class="fc-node-k">03</span><b>Configuration</b><span>paths, steps, sampling, mode<br><code>MonteCarloConfiguration</code></span></li>
<li><span class="fc-node-k">04</span><b>Paths</b><span><code>simulate</code>: N × (steps + 1)</span></li>
<li><span class="fc-node-k">05</span><b>Payoff</b><span>discounted, path by path</span></li>
<li class="fc-out"><span class="fc-node-k">06</span><b>Metrics</b><span><code>pricer</code>, <code>variance</code>, <code>confinter</code>, <code>delta</code>; AD</span></li>
</ol>
<div class="fc-backflow">FinancialFFT.jl, VibratoMonteCarlo.jl and FinancialCalibrator.jl add their own <span>methods</span> to the same <span>pricer</span> <em>← one interface for simulation, Fourier and calibration</em></div>
</div>
```

## Results

```@raw html
<p class="fc-section-lede">Simulated with the package when this site is built.</p>
<div class="fc-tabs">
<input type="radio" name="fc-results" id="fc-results-1" checked><label for="fc-results-1">Paths</label>
<input type="radio" name="fc-results" id="fc-results-2"><label for="fc-results-2">Convergence</label>
<input type="radio" name="fc-results" id="fc-results-3"><label for="fc-results-3">Distributions</label>
<div class="fc-tab-panels">
<figure class="fc-tab-panel"><img src="assets/figures/paths-dark.svg" alt="Simulated paths of three models"><img src="assets/figures/paths-light.svg" alt="Simulated paths of three models"><figcaption><b>Paths.</b> Twenty-four paths of Heston, Kou and Variance Gamma over one year.</figcaption></figure>
<figure class="fc-tab-panel"><img src="assets/figures/convergence-dark.svg" alt="Monte Carlo price and confidence interval against the formula"><img src="assets/figures/convergence-light.svg" alt="Monte Carlo price and confidence interval against the formula"><figcaption><b>Convergence.</b> A Black–Scholes call against the number of paths, with its 95% confidence interval from <code>confinter</code>.</figcaption></figure>
<figure class="fc-tab-panel"><img src="assets/figures/distributions-dark.svg" alt="Density of the log-return at maturity for three models"><img src="assets/figures/distributions-light.svg" alt="Density of the log-return at maturity for three models"><figcaption><b>Distributions.</b> The density of the log-return at one year, from 20 000 paths of each model.</figcaption></figure>
</div>
</div>
```

## Keep in mind

!!! note "Keep in mind"
    First order automatic differentiation works with every payoff, including those for which the
    derivative of the Monte Carlo estimator is meaningless (e.g. binary and barrier options).
    Second order automatic differentiation works too, but the results are seldom useful.
    Multi-currency is not supported.

## Quick start

```@raw html
<p class="fc-section-lede">Five steps, every output computed when this site was built.</p>
```

```@raw html
<div class="fc-step"><span class="fc-step-n">1</span><div><h3>Install the package</h3><p>From the General registry.</p></div></div>
```

```julia
] add FinancialMonteCarlo
```

```@raw html
<div class="fc-step"><span class="fc-step-n">2</span><div><h3>Describe the market and the model</h3><p>σ = 20%, S₀ = 100, a dividend yield of 1%, a flat rate of 2%.</p></div></div>
```

```@example quickstart
using FinancialMonteCarlo

model = BlackScholesProcess(0.2, Underlying(100.0, 0.01))
rfCurve = ZeroRate(0.02)
nothing # hide
```

```@raw html
<div class="fc-step"><span class="fc-step-n">3</span><div><h3>Configure the simulation</h3><p>10 000 paths of 30 steps.</p></div></div>
```

```@example quickstart
mcConfig = MonteCarloConfiguration(10_000, 30)
nothing # hide
```

```@raw html
<div class="fc-step"><span class="fc-step-n">4</span><div><h3>Price an option</h3><p>And its 95% confidence interval.</p></div></div>
```

```@example quickstart
option = EuropeanOption(1.0, 100.0)
(price = pricer(model, rfCurve, mcConfig, option), interval = confinter(model, rfCurve, mcConfig, option, 0.95))
```

```@raw html
<div class="fc-step"><span class="fc-step-n">5</span><div><h3>Price a position</h3><p>Payoffs on the same underlying combine into one contract.</p></div></div>
```

```@example quickstart
position = EuropeanOption(1.0, 100.0) - 2 * BinaryEuropeanOption(1.0, 110.0)
pricer(model, rfCurve, mcConfig, position)
```

```@raw html
<div class="fc-next">
<a href="tutorials/getting_started/"><b>Getting started</b><span>A first tour, with plots.</span></a>
<a href="manual/overview/"><b>Manual</b><span>Every building block in detail.</span></a>
<a href="extending/"><b>Extending the package</b><span>Your own models and payoffs.</span></a>
<a href="api/"><b>API reference</b><span>Every documented name.</span></a>
</div>
```

## The ecosystem

[FinancialFFT.jl](https://github.com/rcalxrc08/FinancialFFT.jl) (Fourier methods),
[VibratoMonteCarlo.jl](https://github.com/rcalxrc08/VibratoMonteCarlo.jl) (smooth Greeks) and
[FinancialCalibrator.jl](https://github.com/rcalxrc08/FinancialCalibrator.jl) (calibration) extend the `pricer` of this
package with methods of their own, on the same models and payoffs. The test suite runs with `] test FinancialMonteCarlo`.
More projects at [rcalxrc08.github.io](https://rcalxrc08.github.io/).
