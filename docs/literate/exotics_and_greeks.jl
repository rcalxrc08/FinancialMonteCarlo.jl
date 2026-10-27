# # Exotic options and greeks
#
# This tutorial goes through the path dependent payoffs and the ways of computing sensitivities.

using FinancialMonteCarlo, DualNumbers

S0, K, r, T, d = 100.0, 100.0, 0.02, 1.0, 0.01
rfCurve = ZeroRate(r)
model = BlackScholesProcess(0.2, Underlying(S0, d))
mcConfig = MonteCarloConfiguration(10_000, 50)

# ## Path dependent payoffs
# Barrier, Asian, American and Bermudan options are priced exactly like European ones.
# American and Bermudan options use the Longstaff-Schwartz regression.
exotics = [
    BarrierOptionDownOut(T, K, 90.0),
    BarrierOptionUpIn(T, K, 120.0),
    DoubleBarrierOption(T, K, 80.0, 130.0),
    AsianFixedStrikeOption(T, K),
    AsianFloatingStrikeOption(T),
    AmericanOption(T, K, false),
    BermudanOption(collect(0.25:0.25:T), K, false),
    BinaryEuropeanOption(T, K),
]
pricer(model, rfCurve, mcConfig, exotics)

# ## Smiles
# A smile is a set of European (or binary) options with the same maturity and different strikes,
# all priced on the same paths.
smile = FinancialMonteCarlo.EuropeanOptionSmile(T, [80.0, 90.0, 100.0, 110.0, 120.0], [false, false, true, true, true])
pricer(model, rfCurve, mcConfig, smile)

# ## Bump and reprice
# `delta` and `rho` compute finite differences with respect to the spot and the zero rate,
# re-using the same seed for both evaluations.
euOption = EuropeanOption(T, K)
delta(model, rfCurve, mcConfig, euOption)
#-
rho(model, rfCurve, mcConfig, euOption)

# ## Automatic differentiation
# Any parameter can be a dual number: the derivative propagates through simulation and payoff.
# Here the dual part of the spot gives the delta,
model_dual_spot = BlackScholesProcess(0.2, Underlying(dual(S0, 1.0), d))
pricer(model_dual_spot, rfCurve, mcConfig, euOption).epsilon

# and the dual part of the volatility gives the vega.
model_dual_vol = BlackScholesProcess(dual(0.2, 1.0), Underlying(S0, d))
pricer(model_dual_vol, rfCurve, mcConfig, euOption).epsilon

# ForwardDiff.jl and ReverseDiff.jl gradients work the same way, see
# [Automatic differentiation](@ref autodiff).
