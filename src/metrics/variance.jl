"""
General Interface for Computation of variance interval of price

		variance_=variance(mcProcess,rfCurve,mcConfig,abstractPayoff)
	
Where:\n
		mcProcess          = Process to be simulated.
		rfCurve  = Zero Rate Data.
		mcConfig = Basic properties of MonteCarlo simulation
		abstractPayoff = Payoff(s) to be priced

		variance_     = variance of the payoff of the derivative

"""
function variance(mcProcess::BaseProcess, rfCurve::AbstractZeroRateCurve, mcConfig::MonteCarloConfiguration, abstractPayoff::AbstractPayoff)
    set_seed(mcConfig)
    T = maturity(abstractPayoff)
    S = simulate(mcProcess, rfCurve, mcConfig, T)
    Payoff = payoff(S, abstractPayoff, rfCurve, mcConfig)
    variance_ = var(Payoff)
    return variance_
end

function variance(mcProcess::BaseProcess, rfCurve::AbstractZeroRateCurve, mcConfig::MonteCarloConfiguration, abstractPayoffs::Array{abstractPayoff_}) where {abstractPayoff_ <: AbstractPayoff}
    set_seed(mcConfig)
    maxT = maximum([maturity(abstractPayoff) for abstractPayoff in abstractPayoffs])
    S = simulate(mcProcess, rfCurve, mcConfig, maxT)
    variance_ = [var(payoff(S, abstractPayoff, rfCurve, mcConfig, maxT)) for abstractPayoff in abstractPayoffs]

    return variance_
end
