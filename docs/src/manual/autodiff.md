```@raw html
<div class="fc-page-head">
<p class="fc-kicker"><span>Manual 10</span> Sensitivities through the simulation</p>
```

# [Automatic differentiation](@id autodiff)

All structs of the package are parametric in the type of their numeric fields, and the numeric
type of the simulated paths is inferred from all the inputs. Any number type implementing the
usual arithmetic can therefore be used for any parameter, and its derivative information
propagates through the simulation and the payoff.

```@raw html
<div class="fc-meta"><span>ForwardDiff</span><span>ReverseDiff</span><span>Mooncake</span><span>Enzyme</span></div>
<nav class="fc-toc">
<a href="#Dual-numbers"><b>01</b>Dual numbers</a>
<a href="#ForwardDiff.jl-and-ReverseDiff.jl"><b>02</b>ForwardDiff.jl and ReverseDiff.jl</a>
<a href="#Mooncake.jl"><b>03</b>Mooncake.jl</a>
<a href="#enzyme"><b>04</b>Enzyme.jl</a>
<a href="#Choosing-a-backend"><b>05</b>Choosing a backend</a>
</nav>
</div>
```

## Dual numbers

```julia
using FinancialMonteCarlo, DualNumbers

model = BlackScholesProcess(0.2)
price = pricer(model, Underlying(dual(100.0, 1.0), 0.01), ZeroRate(0.02), MonteCarloConfiguration(10_000, 30), EuropeanOption(1.0, 100.0))
realpart(price)   # price
price.epsilon     # delta
```

The same works with a dual volatility (vega), rate (rho), maturity (theta) or strike.
[HyperDualNumbers.jl](https://github.com/JuliaDiff/HyperDualNumbers.jl) and
[TaylorSeries.jl](https://github.com/JuliaDiff/TaylorSeries.jl) give higher order derivatives.

## ForwardDiff.jl and ReverseDiff.jl

Write the price as a function of a vector of parameters and take its gradient:

```julia
using ForwardDiff
mcConfig = MonteCarloConfiguration(10_000, 30)
f(x) = pricer(BlackScholesProcess(x[1]), Underlying(x[2], 0.01), ZeroRate(x[3]), mcConfig, EuropeanOption(1.0, 100.0))
ForwardDiff.gradient(f, [0.2, 100.0, 0.02])   # vega, delta, rho
```

`ReverseDiff.gradient` works the same way, but much more slowly (see below). Internally,
validations and other non differentiable bookkeeping are functions marked `ChainRulesCore.@non_differentiable`.

## Mooncake.jl

[Mooncake.jl](https://github.com/chalk-lab/Mooncake.jl) is a reverse mode AD that supports
mutation, so it differentiates directly through the in-place `simulate!` of the models. Seeding
the default `MersenneTwister` ends in a `ccall` into libdSFMT, so the package extension marks
`set_seed!` as having zero derivative; any random number generator works.

```julia
using FinancialMonteCarlo, Mooncake

mcConfig = MonteCarloConfiguration(10_000, 30)

# x = [S0, σ, λ, p, λ₊, λ₋, r, d, T]
f(x) = pricer(KouProcess(x[2], x[3], x[4], x[5], x[6]), Underlying(x[1], x[8]), ZeroRate(x[7]), mcConfig, EuropeanOption(x[9], 100.0))
x0 = [100.0, 0.2, 5.0, 0.3, 30.0, 20.0, 0.02, 0.01, 1.0]

cache = Mooncake.prepare_gradient_cache(f, x0)   # builds the rule, slow the first time
price, (_, grad) = Mooncake.value_and_gradient!!(cache, f, x0)
```

The gradient agrees with `ForwardDiff.gradient(f, x0)` to machine precision. Building the rule
takes a while (about a minute and a half for the call above), after which each gradient costs a
few times a plain pricing; reverse mode pays off when there are many parameters.

## [Enzyme.jl](@id enzyme)

[Enzyme.jl](https://github.com/EnzymeAD/Enzyme.jl) differentiates the compiled code, in forward
and reverse mode, through the in-place simulation as well. With the same `f` and `x0` as above:

```julia
using FinancialMonteCarlo, Enzyme

grad = Enzyme.gradient(Enzyme.Reverse, Enzyme.Const(f), x0)[1]
grad_forward = Enzyme.gradient(Enzyme.Forward, Enzyme.Const(f), x0)[1]   # a tuple, one entry per input
```

Loading Enzyme loads the package extension on EnzymeCore, which marks the seeding of the random
number generator and the scalar draws (jump times and sizes) as inactive: Enzyme cannot analyse the
scalar draws of the `MersenneTwister`, and they carry no derivative. Both modes agree with
ForwardDiff to machine precision for the Black-Scholes, Heston, log-normal mixture, Kou, Merton
and normal inverse Gaussian models, with European, Asian, barrier, American and Bermudan payoffs,
with and without antithetic variates. The variance gamma model does not work: the gamma sampler of
Distributions.jl draws scalar exponentials, which Enzyme cannot analyse either. Positions work too,
also when they are built inside the differentiated function (their legs are two plain vectors);
on Julia 1.10 that case needs runtime activity, `Enzyme.set_runtime_activity(Enzyme.Reverse)`.
Models on a dividend curve need no runtime activity: the nodes of a curve are static vectors, so the
model holds no constant heap memory next to its active parameters (#79).

Zygote is not meant for the simulations, but the model constructors, `set_parameters!`, the curves and
their `integral` are differentiable with it (a calibration loss differentiates through them): `integral`
of a curve and `set_parameters!` have reverse rules (ChainRulesCore), which bring their Zygote gradients
from about 0.5 ms to 4-20 μs (#80).
The first call compiles the derivative, which takes about a minute.

## Choosing a backend

Gradient of a European call with respect to `[S0, σ, r, d]` (for the mixture, the first volatility), 2 000 paths ×
10 steps, minimum time; the first row is the price alone, for comparison. Mooncake's forward mode gives one directional
derivative per call, so a gradient of four parameters costs four times its row:

| Backend | Black-Scholes | Heston | Kou | LogNormal mixture |
|:--|--:|--:|--:|--:|
| price only, `pricer` | 0.27 ms | 0.47 ms | 0.62 ms | 0.87 ms |
| ForwardDiff, `gradient!` (configuration prepared) | 0.75 ms | 1.0 ms | 0.90 ms | 2.3 ms |
| Enzyme, forward (`gradient`, all directions at once) | 0.98 ms | 1.4 ms | 1.3 ms | 2.2 ms |
| Enzyme, reverse | 0.69 ms | 1.0 ms | 2.1 ms | 2.5 ms |
| Mooncake, reverse (`value_and_gradient!!`, rule prepared) | 0.66 ms | 1.2 ms | 2.1 ms | 2.3 ms |
| Mooncake, forward, one direction (`value_and_derivative!!`) | 0.35 ms | 0.65 ms | 1.1 ms | 1.8 ms |

Measured on an Intel Core i7-8550U with Julia 1.12.7, ForwardDiff 1.4.8, Enzyme 0.13.219 and Mooncake 0.5.70.
Before the rules of the kernels below, reverse mode cost, for the same cases, 0.91 / 3.7 / 7.4 / 5.1 ms with Enzyme and
1.5 / 1.5 / 2.4 / 8.9 ms with Mooncake; an Asian payoff made Enzyme allocate once per path.

Both reverse modes run within a few times the price because the loops the models spend their time in are functions of
their own with hand-written derivatives: a step of the Brownian motion, of Heston and of the mixture, the exponential of
the log paths, the carry factors of Bachelier and of the shifted models, and the jumps of Kou and Merton, which are drawn
first as random numbers and grid points, and sized from the parameters after. The rules work on whole arrays of
`Float64` (or `Float32`); without them Mooncake records every element of every broadcast and Enzyme caches every value
the paths overwrite. Dual numbers and the GPU use the generic code.

- ForwardDiff carries all the partials through one simulation: it is the simplest choice up to a few dozen
  parameters.
- Enzyme and Mooncake in reverse mode cost the same whatever the number of parameters, the choice for many of them
  (the calibration of a large model). Mooncake also works with the variance gamma model, whose gamma sampler Enzyme
  cannot analyse.
- Mooncake's forward mode, one direction at a time, is the cheapest single sensitivity: about the price itself.
- ReverseDiff records one tape entry per number of every path, two orders of magnitude slower than ForwardDiff (65 ms
  for the Black-Scholes case, 12 ms with a compiled tape, Julia 1.13.3): a compiled tape
  (`ReverseDiff.compile(ReverseDiff.GradientTape(f, x0))`) replays the operations recorded at `x0`, which is correct here
  because the seed is reset at every call, but not if the number of paths, the seed or the payoff changes, nor when a
  parameter changes the control flow (the number of jumps, the exercise decisions of an American option).
- TaylorSeries allocates a vector of coefficients for every number of every path (5 ms for a vega of order 2).

With `MultiThreading` and `MultiProcess` only the backends based on number types (ForwardDiff,
DualNumbers, HyperDualNumbers, TaylorSeries) are supported.

!!! warning
    The derivative of the Monte Carlo estimator is the estimator of the derivative only for
    payoffs that are (almost everywhere) differentiable along the paths. Automatic differentiation
    runs for every payoff, but for discontinuous ones (binary, barrier) the result is not the greek.
    Second order derivatives run as well, but are seldom meaningful.

    The same applies to parameters that change the number of jumps in a path. For the jump
    diffusions (Kou, Merton) the jump times are drawn as `randexp() / (λ * T)`, so the derivatives
    with respect to `T` and `λ` miss the contribution of jumps crossing the maturity; the other
    sensitivities are unaffected. With `StandardMC` the direction of a Kou jump is drawn against `p`,
    so the derivative with respect to `p` misses its effect on the jumps; with `AntitheticMC` and
    `SobolMC` the sizes are quantiles of uniforms, continuous in `p`, and the derivative is complete.
