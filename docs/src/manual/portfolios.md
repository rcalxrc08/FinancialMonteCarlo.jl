# Portfolios and multivariate

## Positions: algebra of contracts

Contracts on the same underlying form a vector space over the real numbers. Adding, subtracting
or scaling payoffs returns a [`FinancialMonteCarlo.Position`](@ref): the payoffs with their
weights, in the order they were added.

```julia
straddle = EuropeanOption(T, K) + EuropeanOption(T, K, false)
spread = EuropeanOption(T, 95.0) - EuropeanOption(T, 105.0)
position = 2.0 * straddle - spread / 2
pricer(model, rfCurve, mcConfig, position)
```

- The same contract appears once and its weights add up. Two payoffs are the same contract when
  they have the same type and equal parameters, whatever their number types:
  `EuropeanOption(T, 100) + EuropeanOption(T, 100.0)` is one leg of weight 2.
- A leg whose weight is exactly zero is removed: `c - c` is the empty position, and a netted leg no
  longer sets how far the paths are simulated. A weight that is a dual number with value zero is
  kept, since it can carry a derivative.
- The weights are promoted to one type, so a position of payoffs of one type (spreads, butterflies,
  strips) is concretely typed and its price is inferred.
- `position[payoff]` is the weight of `payoff` (zero if absent), and a position iterates as
  `payoff => weight` pairs. `FinancialMonteCarlo.Position(payoff => weight, ...)` builds one
  directly, and `sum(positions)` adds several.

A position is priced on a single simulation, up to the largest maturity of its legs, path by path:
its discounted payoff along each path is the weighted sum of those of its legs. So `variance` and
`confinter` of a position are those of the position as a whole, with the correlation between its
legs, and `delta` and `rho` bump and reprice it like a single payoff.

## Market data sets and portfolios

To price contracts on several underlyings, give each underlying a name and link it with `→`
(typed `\to<TAB>`):

- a [`MarketDataSet`](@ref) maps names to models: `"ENI" → model`;
- a [`FinancialMonteCarlo.Portfolio`](@ref) maps names to positions: `"ENI" → payoff` or
  `"ENI" → position`.

Both can be combined with `+`; adding portfolios adds the positions on the same name.

```julia
market = ("ENI" → model_eni) + ("AAPL" → model_aapl)
portfolio = ("ENI" → straddle) + ("AAPL" → 3.0 * AsianFixedStrikeOption(T, K))
pricer(market, rfCurve, mcConfig, portfolio)
delta(market, rfCurve, mcConfig, portfolio, "ENI")
```

Each underlying is simulated once. The models are simulated one after the other, in the order of
the market data set, from one stream of random numbers seeded once, so the underlyings are
independent (with `SobolMC`, which would give every model the same points, a portfolio on several
models is an error). Names in the market data set that have no position are not simulated.

!!! warning
    Positions on names that are missing from the market data set are not priced, i.e. they
    contribute zero to the price.

## Multivariate models

A multivariate model is linked to the tuple of the names of its components:

```julia
ρ = [1.0 0.3; 0.3 1.0]
market = ("ENI", "AAPL") → GaussianCopulaNVariateProcess(ρ, model_eni, model_aapl)
```

A portfolio on such a market can contain

- **basket payoffs** (e.g. [`EuropeanOptionND`](@ref)) on the whole basket: `("ENI", "AAPL") → EuropeanOptionND(T, K)`;
- **single name payoffs** on one of the components: `"ENI" → EuropeanOption(T, K)`, priced on
  the marginal of the joint simulation.

Names joined by `_` (`"ENI_AAPL"`) still work, with a deprecation warning, and become the tuple of
the names.

See [Portfolios](@ref api-portfolios) in the API reference.
