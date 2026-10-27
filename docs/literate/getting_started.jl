# # Getting started
#
# This tutorial prices a European call under Black-Scholes and shows the handful of objects every
# computation in FinancialMonteCarlo.jl is made of.
#
# ## Installation
# The package is registered, so from the Julia REPL:
# ```julia
# ] add FinancialMonteCarlo
# ```

using FinancialMonteCarlo

# ## The four ingredients
# Every metric takes the same four arguments: a **model**, a **zero rate curve**, a
# **Monte Carlo configuration** and a **payoff** (or a collection of them).
#
# The market data: spot price, continuous dividend yield and risk free rate.
S0 = 100.0
d = 0.01
r = 0.02
underlying = Underlying(S0, d)
rfCurve = ZeroRate(r)

# The model of the underlying, here Black-Scholes with 20% volatility:
model = BlackScholesProcess(0.2, underlying)

# The simulation settings: number of paths and number of time steps per path.
mcConfig = MonteCarloConfiguration(10_000, 30)

# The contract: a European call with maturity one year and strike 100.
T = 1.0
K = 100.0
euOption = EuropeanOption(T, K)

# ## Pricing
pricer(model, rfCurve, mcConfig, euOption)

# Prices are reproducible: the configuration carries a seed, and every metric resets the random
# number generator before simulating.
pricer(model, rfCurve, mcConfig, euOption) == pricer(model, rfCurve, mcConfig, euOption)

# Several payoffs can be priced on the same simulated paths by passing a vector:
payoffs = [Forward(T), EuropeanOption(T, K), EuropeanOption(T, K, false), AmericanOption(T, K, false)]
pricer(model, rfCurve, mcConfig, payoffs)

# ## Accuracy of the estimate
# `variance` and `confinter` return the variance of the discounted payoff and a confidence interval
# for the price.
variance(model, rfCurve, mcConfig, euOption)
#-
confinter(model, rfCurve, mcConfig, euOption, 0.99)

# Antithetic variates reduce the variance at no extra cost (the number of paths must be even):
mcAntithetic = MonteCarloConfiguration(10_000, 30, FinancialMonteCarlo.AntitheticMC())
variance(model, rfCurve, mcAntithetic, euOption)

# ## Changing the model
# Models are interchangeable: the same payoffs can be priced under, e.g., Heston or Kou.
heston = HestonProcess(0.2, 0.2, 0.01, 0.01, 0.0, 0.03, underlying)
kou = KouProcess(0.2, 5.0, 0.3, 30.0, 20.0, underlying)
[pricer(m, rfCurve, mcConfig, euOption) for m in (model, heston, kou)]

# ## Looking at the paths
# `simulate` returns the simulated paths as an `Nsim × (Nstep + 1)` matrix, the first column being the spot.
S = simulate(model, rfCurve, mcConfig, T)
size(S)
