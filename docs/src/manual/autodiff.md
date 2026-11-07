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

## Choosing a backend

Gradient of a European call with respect to `[S0, σ, r, d]`, 2 000 paths × 10 steps, minimum time:

| Backend | Black-Scholes | Kou | allocations (Black-Scholes) |
|:--|--:|--:|--:|
| ForwardDiff, `gradient` | 0.53 ms | 0.97 ms | 29 |
| Mooncake, `value_and_gradient!!` (rule prepared) | 2.5 ms | 3.9 ms | 96 |
| ReverseDiff, `gradient` | 68 ms | 81 ms | 576 220 |
| TaylorSeries, vega of order 2 | 4.0 ms | | 120 047 |

Measured on an Intel Core i7-8550U with Julia 1.13. ForwardDiff carries all the partials through one
simulation, so it is the fastest up to a few dozen parameters; Mooncake pays a fixed cost per
gradient and is the choice for many parameters, such as the calibration of a large model.
ReverseDiff records one tape entry per number of every path, which makes it two orders of
magnitude slower than ForwardDiff here; prefer Mooncake for reverse mode. With `MultiThreading`
and `MultiProcess` only the forward backends are supported.

!!! warning
    The derivative of the Monte Carlo estimator is the estimator of the derivative only for
    payoffs that are (almost everywhere) differentiable along the paths. Automatic differentiation
    runs for every payoff, but for discontinuous ones (binary, barrier) the result is not the greek.
    Second order derivatives run as well, but are seldom meaningful.

    The same applies to parameters that change the number of jumps in a path. For the jump
    diffusions (Kou, Merton) the jump times are drawn as `randexp() / (λ * T)`, so the derivatives
    with respect to `T` and `λ` miss the contribution of jumps crossing the maturity; the other
    sensitivities are unaffected.
