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

`ReverseDiff.gradient` works the same way. Internally, assertions and other non differentiable
bookkeeping are wrapped in `ChainRulesCore.@ignore_derivatives`.

## Mooncake.jl

[Mooncake.jl](https://github.com/chalk-lab/Mooncake.jl) is a reverse mode AD that supports
mutation, so it differentiates directly through the in-place `simulate!` of the models. The only
requirement is a pure Julia random number generator: the default `MersenneTwister` seeds and
refills itself through `ccall`s into libdSFMT, for which Mooncake has no rule, while `Xoshiro`
works out of the box.

```julia
using FinancialMonteCarlo, Mooncake, Random

mcConfig = MonteCarloConfiguration(10_000, 30, FinancialMonteCarlo.StandardMC(), FinancialMonteCarlo.SerialMode(0, Xoshiro()))

# x = [S0, σ, λ, p, λ₊, λ₋, r, d, T]
f(x) = pricer(KouProcess(x[2], x[3], x[4], x[5], x[6], Underlying(x[1], x[8])), ZeroRate(x[7]), mcConfig, EuropeanOption(x[9], 100.0))
x0 = [100.0, 0.2, 5.0, 0.3, 30.0, 20.0, 0.02, 0.01, 1.0]

cache = Mooncake.prepare_gradient_cache(f, x0)   # builds the rule, slow the first time
price, (_, grad) = Mooncake.value_and_gradient!!(cache, f, x0)
```

The gradient agrees with `ForwardDiff.gradient(f, x0)` to machine precision. Building the rule
takes a while (about a minute and a half for the call above), after which each gradient costs a
few times a plain pricing; reverse mode pays off when there are many parameters.

!!! warning
    The derivative of the Monte Carlo estimator is the estimator of the derivative only for
    payoffs that are (almost everywhere) differentiable along the paths. Automatic differentiation
    runs for every payoff, but for discontinuous ones (binary, barrier) the result is not the greek.
    Second order derivatives run as well, but are seldom meaningful.

    The same applies to parameters that change the number of jumps in a path. For the jump
    diffusions (Kou, Merton) the jump times are drawn as `randexp() / (λ * T)`, so the derivatives
    with respect to `T` and `λ` miss the contribution of jumps crossing the maturity; the other
    sensitivities are unaffected.
