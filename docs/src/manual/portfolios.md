# Portfolios and multivariate

## Positions: algebra of contracts

Contracts on the same underlying form a vector space over the real numbers. Adding, subtracting
or scaling payoffs returns a `Position`, a dictionary from payoff to weight:

```julia
straddle = EuropeanOption(T, K) + EuropeanOption(T, K, false)
spread = EuropeanOption(T, 95.0) - EuropeanOption(T, 105.0)
position = 2.0 * straddle - spread / 2
pricer(model, rfCurve, mcConfig, position)
```

A position is priced on a single simulation, up to the largest maturity of its payoffs.

## Market data sets and portfolios

To price contracts on several underlyings, give each underlying a name and link it with `→`
(typed `\to<TAB>`):

- a [`MarketDataSet`](@ref) maps names to models: `"ENI" → model`;
- a portfolio maps names to positions: `"ENI" → payoff` or `"ENI" → position`.

Both can be combined with `+`:

```julia
market = ("ENI" → model_eni) + ("AAPL" → model_aapl)
portfolio = ("ENI" → straddle) + ("AAPL" → 3.0 * AsianFixedStrikeOption(T, K))
pricer(market, rfCurve, mcConfig, portfolio)
```

Each underlying is simulated once, and names in the market data set that do not appear in the
portfolio are not simulated at all. Single name keys cannot contain `_`.

!!! warning
    Contracts on names that are missing from the market data set are silently skipped, i.e. they
    contribute zero to the price.

## Multivariate models

The name of a multivariate model is the list of the names of its components joined by `_`:

```julia
ρ = [1.0 0.3; 0.3 1.0]
market = "ENI_AAPL" → GaussianCopulaNVariateProcess(ρ, model_eni, model_aapl)
```

A portfolio on such a market can contain

- **basket payoffs** (e.g. [`EuropeanOptionND`](@ref)) on the whole basket: `"ENI_AAPL" → EuropeanOptionND(T, K)`;
- **single name payoffs** on one of the components: `"ENI" → EuropeanOption(T, K)`, priced on
  the marginal of the joint simulation.

See [Portfolios](@ref api-portfolios) in the API reference.
