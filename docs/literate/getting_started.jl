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

# The interval shrinks like one over the square root of the number of paths. Under Black-Scholes the
# exact price is known, so the convergence can be checked (a European payoff needs one time step):
using Plots, Distributions
default(size = (720, 400), linewidth = 2, framestyle = :axes, grid = :y, gridalpha = 0.25, titlefontsize = 11, guidefontsize = 10, legendfontsize = 9)

σ = 0.2
d1 = (log(S0 / K) + (r - d + σ^2 / 2) * T) / (σ * sqrt(T))
exact = S0 * exp(-d * T) * cdf(Normal(), d1) - K * exp(-r * T) * cdf(Normal(), d1 - σ * sqrt(T))
paths = [500, 1_000, 2_000, 5_000, 10_000, 20_000, 50_000, 100_000]
intervals = [confinter(model, rfCurve, MonteCarloConfiguration(n, 1), euOption, 0.99) for n in paths]
centres = [(lower + upper) / 2 for (lower, upper) in intervals]
plot(paths, centres; ribbon = (centres .- first.(intervals), last.(intervals) .- centres), fillalpha = 0.2, color = "#2a78d6", marker = :circle, markersize = 4, markerstrokewidth = 0, label = "Monte Carlo, 99% interval", xscale = :log10, xticks = ([1_000, 10_000, 100_000], ["1k", "10k", "100k"]), xlabel = "number of paths", ylabel = "price", title = "European call: Monte Carlo estimate against the exact price")
hline!([exact]; color = "#eb6834", linestyle = :dash, label = "exact Black-Scholes price")

# ## Changing the model
# Models are interchangeable: the same payoffs can be priced under, e.g., Heston or Kou.
heston = HestonProcess(0.2, 0.2, 0.01, 0.01, 0.0, 0.03, underlying)
kou = KouProcess(0.2, 5.0, 0.3, 30.0, 20.0, underlying)
[pricer(m, rfCurve, mcConfig, euOption) for m in (model, heston, kou)]

# ## Looking at the paths
# `simulate` returns the simulated paths as an `Nsim × (Nstep + 1)` matrix, the first column being the spot.
S = simulate(model, rfCurve, mcConfig, T)
size(S)

# The first paths, with the strike of the call:
times = range(0, T, length = size(S, 2))
plot(times, S[1:30, :]'; color = "#2a78d6", alpha = 0.35, linewidth = 1, label = ["paths" fill("", 1, 29)], legend = :topleft, xlabel = "time (years)", ylabel = "underlying", title = "30 simulated Black-Scholes paths")
hline!([K]; color = :gray40, linestyle = :dash, label = "strike")
