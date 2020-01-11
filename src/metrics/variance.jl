"""
General Interface for Computation of variance interval of price

		variance_=variance(mcProcess,rfCurve,mcBaseData,payoff_,)
	
Where:\n
		mcProcess          = Process to be simulated.
		rfCurve  = Datas of the Spot.
		mcBaseData = Basic properties of MonteCarlo simulation
		payoff_ = Payoff(s) to be priced
		

		variance_     = variance of the payoff of the derivative

"""
function variance(mcProcess::BaseProcess,rfCurve::AbstractZeroRateCurve,mcConfig::MonteCarloConfiguration,abstractPayoff::AbstractPayoff)
	set_seed(mcConfig)
	T=maturity(abstractPayoff);
	S=simulate(mcProcess,rfCurve,mcConfig,T)
	Payoff=payoff(S,abstractPayoff,rfCurve);
	variance_=var(Payoff);
	return variance_;
end


function variance(mcProcess::BaseProcess,rfCurve::AbstractZeroRateCurve,mcConfig::MonteCarloConfiguration,abstractPayoffs::Array{abstractPayoff_}) where {abstractPayoff_ <: AbstractPayoff}
	set_seed(mcConfig)
	maxT=maximum([maturity(abstractPayoff) for abstractPayoff in abstractPayoffs])
	S=simulate(mcProcess,rfCurve,mcConfig,maxT)
	variance_=[var(payoff(S,abstractPayoff,rfCurve,maxT)) for abstractPayoff in abstractPayoffs  ]
	
	return variance_;
end
