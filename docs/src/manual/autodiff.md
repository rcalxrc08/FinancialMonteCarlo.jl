# [Automatic differentiation](@id autodiff)

All structs of the package are parametric in the type of their numeric fields, and the numeric
type of the simulated paths is inferred from all the inputs. Any number type implementing the
usual arithmetic can therefore be used for any parameter, and its derivative information
propagates through the simulation and the payoff.

## Dual numbers

```julia
using FinancialMonteCarlo, DualNumbers

model = BlackScholesProcess(0.2, Underlying(dual(100.0, 1.0), 0.01))
price = pricer(model, ZeroRate(0.02), MonteCarloConfiguration(10_000, 30), EuropeanOption(1.0, 100.0))
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
f(x) = pricer(BlackScholesProcess(x[1], Underlying(x[2], 0.01)), ZeroRate(x[3]), mcConfig, EuropeanOption(1.0, 100.0))
ForwardDiff.gradient(f, [0.2, 100.0, 0.02])   # vega, delta, rho
```

`ReverseDiff.gradient` works the same way, but much more slowly (see below). Internally,
assertions and other non differentiable bookkeeping are wrapped in `ChainRulesCore.@ignore_derivatives`.

## Mooncake.jl

[Mooncake.jl](https://github.com/chalk-lab/Mooncake.jl) is a reverse mode AD that supports
mutation, so it differentiates directly through the in-place `simulate!` of the models. Seeding
the default `MersenneTwister` ends in a `ccall` into libdSFMT, so the package extension marks
`set_seed!` as having zero derivative; any random number generator works.

```julia
using FinancialMonteCarlo, Mooncake

mcConfig = MonteCarloConfiguration(10_000, 30)

# x = [S0, σ, λ, p, λ₊, λ₋, r, d, T]
f(x) = pricer(KouProcess(x[2], x[3], x[4], x[5], x[6], Underlying(x[1], x[8])), ZeroRate(x[7]), mcConfig, EuropeanOption(x[9], 100.0))
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
also when they are built inside the differentiated function (their legs are two plain vectors).
The first call compiles the derivative, which takes about a minute.

## Choosing a backend

Gradient of a European call with respect to `[S0, σ, r, d]`, 2 000 paths × 10 steps, minimum time
over three runs; the first row is the price alone, for comparison:

| Backend | Black-Scholes | Kou | allocations (Black-Scholes) | allocations (Kou) |
|:--|--:|--:|--:|--:|
| price only, `pricer` | 0.27 ms | 0.97 ms | 21 | 21 |
| ForwardDiff, `gradient!` (configuration prepared) | 0.49 ms | 1.1 ms | 21 | 21 |
| Enzyme, forward | 0.63 ms | 1.3 ms | 77 | 77 |
| Enzyme, reverse | 0.90 ms | 7.7 ms | 38 | 5 036 |
| Mooncake, `value_and_gradient!!` (rule prepared) | 2.3 ms | 4.5 ms | 90 | 90 |
| ReverseDiff, `gradient!` with a compiled tape | 12 ms | 19 ms | 0 | 0 |
| ReverseDiff, `gradient` | 65 ms | 104 ms | 596 198 | 928 109 |
| TaylorSeries, vega of order 2 | 5.0 ms | | 120 041 | |

Measured on an Intel Core i7-8550U with Julia 1.13.3, ForwardDiff 1.4.7, Enzyme 0.13.219,
Mooncake 0.5.69, ReverseDiff 1.18.7 and TaylorSeries 0.23.3. The price allocates the paths, the
normal draws of a step and the discounted payoffs, plus the eleven small allocations of
`Random.seed!` for the `MersenneTwister`.

- ForwardDiff carries all the partials through one simulation, at less than twice the price
  here: it is the fastest up to a few dozen parameters.
- Enzyme is the fastest reverse mode: in reverse its cost does not grow with the number of
  parameters, so it is the choice for many of them, such as the calibration of a large model. For
  the jump models it caches the jumps of every path (a loop whose length is random), which makes
  it slower than Mooncake there.
- Mooncake pays a fixed cost per gradient and works with every model, variance gamma included.
- ReverseDiff records one tape entry per number of every path, which makes it two orders of
  magnitude slower than ForwardDiff. A compiled tape (`ReverseDiff.compile(ReverseDiff.GradientTape(f, x0))`)
  is five times faster and allocates nothing; it replays the operations recorded at `x0`, which
  is correct here because the seed is reset at every call, but not if the number of paths, the
  seed or the payoff changes, nor when a parameter changes the control flow (the number of jumps,
  the exercise decisions of an American option).
- TaylorSeries allocates a vector of coefficients for every number of every path.

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
