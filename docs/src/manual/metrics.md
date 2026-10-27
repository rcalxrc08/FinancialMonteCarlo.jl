# Metrics

Metrics are functions of a model, a zero rate curve, a Monte Carlo configuration and a payoff.
They all accept a single payoff or a vector of payoffs (priced on the same paths), and most of
them also positions and portfolios.

| Metric | Returns |
|:--|:--|
| [`pricer`](@ref) | Monte Carlo estimate of the price |
| [`variance`](@ref) | variance of the discounted payoff |
| [`confinter`](@ref) | confidence interval of the price at a given level (default 99%) |
| [`delta`](@ref) | sensitivity to the spot price, by bump and reprice |
| [`rho`](@ref) | sensitivity to the zero rate, by bump and reprice |
| [`FinancialMonteCarlo.distribution`](@ref) | simulated values of the underlying at a given time |

```julia
pricer(model, rfCurve, mcConfig, payoff)
pricer(model, rfCurve, mcConfig, [payoff_1, payoff_2])     # vector of prices
pricer(model, rfCurve, mcConfig, 2.0 * payoff_1 - payoff_2) # price of a position
variance(model, rfCurve, mcConfig, payoff)
confinter(model, rfCurve, mcConfig, payoff, 0.95)
delta(model, rfCurve, mcConfig, payoff, 1e-4)               # spot bump
rho(model, rfCurve, mcConfig, payoff, 1e-4)                 # rate bump
```

The bump and reprice greeks use the same seed for the base and the bumped prices, so the Monte
Carlo noise largely cancels. For exact pathwise derivatives, use dual numbers instead, see
[Automatic differentiation](@ref autodiff).

For a market data set and a portfolio, [`delta`](@ref) takes the name of the underlying to bump:

```julia
delta(market, rfCurve, mcConfig, portfolio, "ENI")
```

See [Metrics](@ref api-metrics) in the API reference.
