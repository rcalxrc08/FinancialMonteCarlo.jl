# # Exotic options and greeks
#
# This tutorial goes through the path dependent payoffs and the ways of computing sensitivities.

using FinancialMonteCarlo, DualNumbers, Plots
default(size = (720, 400), linewidth = 2, framestyle = :axes, grid = :y, gridalpha = 0.25, titlefontsize = 11, guidefontsize = 10, legendfontsize = 9)

S0, K, r, T, d = 100.0, 100.0, 0.02, 1.0, 0.01
rfCurve = ZeroRate(r)
model = BlackScholesProcess(0.2)
und = Underlying(S0, d) # spot and dividends, passed next to the model
mcConfig = MonteCarloConfiguration(10_000, 50)

# ## Path dependent payoffs
# Barrier, Asian, American and Bermudan options are priced exactly like European ones.
# American and Bermudan options use the Longstaff-Schwartz regression.
exotics = [BarrierOptionDownOut(T, K, 90.0), BarrierOptionUpIn(T, K, 120.0), DoubleBarrierOption(T, K, 80.0, 130.0), AsianFixedStrikeOption(T, K), AsianFloatingStrikeOption(T), AmericanOption(T, K, false), BermudanOption(collect(0.25:0.25:T), K, false), BinaryEuropeanOption(T, K)]
pricer(model, und, rfCurve, mcConfig, exotics)

# ## Smiles
# A smile is a set of European (or binary) options with the same maturity and different strikes,
# all priced on the same paths.
smile = FinancialMonteCarlo.EuropeanOptionSmile(T, [80.0, 90.0, 100.0, 110.0, 120.0], [false, false, true, true, true])
pricer(model, und, rfCurve, mcConfig, smile)

# A smile prices many strikes on a single simulation, so a whole price curve costs one run:
strikes = collect(70.0:2.0:130.0)
calls = pricer(model, und, rfCurve, mcConfig, FinancialMonteCarlo.EuropeanOptionSmile(T, strikes, fill(true, length(strikes))))
puts = pricer(model, und, rfCurve, mcConfig, FinancialMonteCarlo.EuropeanOptionSmile(T, strikes, fill(false, length(strikes))))
plot(strikes, calls; color = "#2a78d6", label = "calls", legend = :top, xlabel = "strike", ylabel = "price", title = "European options by strike, one simulation per curve")
plot!(strikes, puts; color = "#eb6834", label = "puts")
annotate!([(strikes[3], calls[3], text("calls", 9, "#2a78d6", :left, :bottom)), (strikes[end-5], puts[end-5], text("puts", 9, "#eb6834", :right, :bottom))])

# ## Implied volatility
# Prices of different strikes are easier to compare as Black-Scholes implied volatilities: the volatility that,
# in the Black-Scholes formula, gives the Monte Carlo price. `blsimpv` from
# [FinancialToolbox.jl](https://github.com/rcalxrc08/FinancialToolbox.jl) inverts the formula. Black-Scholes itself
# gives a flat line at its volatility; the other models give a smile or a skew. Each curve is priced with one smile
# of out of the money options (puts below the spot, calls above), which have the smallest Monte Carlo error.
using FinancialToolbox: blsimpv

smile_strikes = collect(70.0:5.0:140.0)
otm_is_call = smile_strikes .>= S0
otm_smile = FinancialMonteCarlo.EuropeanOptionSmile(T, smile_strikes, otm_is_call)
mcSmile = MonteCarloConfiguration(100_000, 100; sampling = FinancialMonteCarlo.AntitheticMC())
smile_models = (("Black-Scholes", model, "#2a78d6"), ("Heston, ρ = -0.7", HestonProcess(0.5, 0.2, 0.0, 2.0, -0.7, 0.04), "#eb6834"), ("Kou", KouProcess(0.2, 5.0, 0.3, 30.0, 20.0), "#1baf7a"), ("Variance Gamma", VarianceGammaProcess(0.2, -0.1, 0.2), "#eda100"))
plot(; legend = :bottomleft, xlabel = "strike", ylabel = "implied volatility", title = "Implied volatility by strike, one year maturity")
for (name, m, colour) in smile_models
    prices = pricer(m, und, rfCurve, mcSmile, otm_smile)
    plot!(smile_strikes, blsimpv.(S0, smile_strikes, r, T, prices, d, otm_is_call); color = colour, marker = :circle, markersize = 3, markerstrokewidth = 0, label = name)
end
plot!()

# ## Bump and reprice
# `delta` and `rho` compute finite differences with respect to the spot and the zero rate,
# re-using the same seed for both evaluations.
euOption = EuropeanOption(T, K)
delta(model, und, rfCurve, mcConfig, euOption)
#-
rho(model, und, rfCurve, mcConfig, euOption)

# The delta of a barrier option can behave very differently from the one of a vanilla option. An up-and-out call
# gains from a rise of the spot only while the barrier is far: closer to it, a rise makes the knock out more
# likely, and the delta falls towards zero while the one of the vanilla call keeps growing. (The barrier is
# checked only at the 50 simulation dates, which softens the fall: with continuous monitoring the delta would
# turn negative next to the barrier.)
barrier_option = BarrierOptionUpOut(T, K, 130.0)
spots = collect(70.0:1.0:129.0)
delta_vanilla = [delta(model, Underlying(S, d), rfCurve, mcConfig, euOption) for S in spots]
delta_barrier = [delta(model, Underlying(S, d), rfCurve, mcConfig, barrier_option) for S in spots]
plot(spots, delta_vanilla; color = "#2a78d6", label = "European call", legend = :topleft, xlabel = "spot", ylabel = "delta", title = "Delta by spot: European call and up-and-out call (barrier 130)")
plot!(spots, delta_barrier; color = "#eb6834", label = "up-and-out call")
vline!([130.0]; color = :gray40, linestyle = :dash, linewidth = 1, label = "barrier")

# ## Automatic differentiation
# Any parameter can be a dual number: the derivative propagates through simulation and payoff.
# Here the dual part of the spot gives the delta,
pricer(model, Underlying(dual(S0, 1.0), d), rfCurve, mcConfig, euOption).epsilon

# and the dual part of the volatility gives the vega.
model_dual_vol = BlackScholesProcess(dual(0.2, 1.0))
pricer(model_dual_vol, und, rfCurve, mcConfig, euOption).epsilon

# ForwardDiff.jl and ReverseDiff.jl gradients work the same way, see
# [Automatic differentiation](@ref autodiff).
