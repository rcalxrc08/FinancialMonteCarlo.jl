# # Portfolios and multivariate models
#
# Contracts can be combined into positions and portfolios, and priced against a set of models,
# one for each underlying.

using FinancialMonteCarlo

S0, K, r, T, d = 100.0, 100.0, 0.02, 1.0, 0.01
rfCurve = ZeroRate(r)
mcConfig = MonteCarloConfiguration(10_000, 30)

# ## Positions
# Payoffs on the same underlying form a vector space: they can be added and scaled.
# The result is a `Position`, i.e. a dictionary from payoff to weight.
straddle = EuropeanOption(T, K) + EuropeanOption(T, K, false)
#-
bull_spread = EuropeanOption(T, 95.0) - EuropeanOption(T, 105.0)
#-
model = BlackScholesProcess(0.2, Underlying(S0, d))
pricer(model, rfCurve, mcConfig, 2.0 * straddle - bull_spread)

# ## Market data sets and portfolios
# A `MarketDataSet` links underlying names to models with the `→` operator (typed `\to<TAB>`),
# a portfolio links underlying names to positions.
market = "ENI" → BlackScholesProcess(0.2, Underlying(S0, d))
market += "AAPL" → KouProcess(0.25, 5.0, 0.3, 30.0, 20.0, Underlying(150.0, 0.0))

portfolio = "ENI" → straddle
portfolio += "AAPL" → 3.0 * AsianFixedStrikeOption(T, 150.0)

pricer(market, rfCurve, mcConfig, portfolio)

# Delta with respect to a single underlying of the market data set:
delta(market, rfCurve, mcConfig, portfolio, "AAPL")

# ## Multivariate models
# Correlated underlyings are modelled with a Gaussian copula on top of single name models.
# The key of a multivariate model joins the names of its components with `_`.
ρ = [1.0 0.5; 0.5 1.0]
basket_model = GaussianCopulaNVariateProcess(ρ, BlackScholesProcess(0.2, Underlying(S0, d)), BlackScholesProcess(0.3, Underlying(S0, d)))
basket_market = "ENI_AAPL" → basket_model

# Basket payoffs such as `EuropeanOptionND` refer to the whole basket, while single name payoffs
# can still refer to a single component.
basket_portfolio = "ENI_AAPL" → EuropeanOptionND(T, K)
basket_portfolio += "ENI" → EuropeanOption(T, K)
pricer(basket_market, rfCurve, mcConfig, basket_portfolio)
