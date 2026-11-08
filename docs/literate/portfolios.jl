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

# A position is priced like a single payoff, so its value can be drawn against the spot, here with its
# payoff at maturity (the value it converges to as the maturity approaches):
using Plots
default(size = (720, 400), linewidth = 2, framestyle = :axes, grid = :y, gridalpha = 0.25, titlefontsize = 11, guidefontsize = 10, legendfontsize = 9)

spots = collect(70.0:2.0:130.0)
value_today(position, S) = pricer(BlackScholesProcess(0.2, Underlying(S, d)), rfCurve, mcConfig, position)
plot(spots, [value_today(straddle, S) for S in spots]; color = "#2a78d6", label = "straddle, today", legend = :top, xlabel = "spot", ylabel = "value", title = "Value of two positions by spot, today and at maturity")
plot!(spots, abs.(spots .- K); color = "#2a78d6", linestyle = :dash, linewidth = 1, label = "straddle, at maturity")
plot!(spots, [value_today(bull_spread, S) for S in spots]; color = "#eb6834", label = "bull spread 95-105, today")
plot!(spots, clamp.(spots .- 95.0, 0.0, 10.0); color = "#eb6834", linestyle = :dash, linewidth = 1, label = "bull spread, at maturity")

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

# `simulate` returns one path matrix per component. The log returns at maturity of the two components are
# correlated as the copula prescribes, while each keeps the distribution of its own model (the second is wider,
# since its volatility is 30%):
S_basket = simulate(basket_model, rfCurve, mcConfig, T)
returns_1 = log.(S_basket[1][:, end] ./ S0)
returns_2 = log.(S_basket[2][:, end] ./ S0)
scatter(returns_1[1:2_000], returns_2[1:2_000]; color = "#2a78d6", markersize = 3, markerstrokewidth = 0, alpha = 0.4, label = "", xlabel = "log return of ENI", ylabel = "log return of AAPL", title = "Log returns at maturity of a Gaussian copula with ρ = 0.5 (2000 paths)")

# Basket payoffs such as `EuropeanOptionND` refer to the whole basket, while single name payoffs
# can still refer to a single component.
basket_portfolio = "ENI_AAPL" → EuropeanOptionND(T, K)
basket_portfolio += "ENI" → EuropeanOption(T, K)
pricer(basket_market, rfCurve, mcConfig, basket_portfolio)
