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

!!! warning
    The derivative of the Monte Carlo estimator is the estimator of the derivative only for
    payoffs that are (almost everywhere) differentiable along the paths. Automatic differentiation
    runs for every payoff, but for discontinuous ones (binary, barrier) the result is not the greek.
    Second order derivatives run as well, but are seldom meaningful.
