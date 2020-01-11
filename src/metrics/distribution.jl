"""
General Interface for Pricing

		Price=pricer(mcProcess,rfCurve,mcBaseData,payoff_)
	
Where:\n
		mcProcess          = Process to be simulated.
		rfCurve  = Datas of the Spot.
		mcBaseData = Basic properties of MonteCarlo simulation
		payoff_ = Payoff(s) to be priced
		

		Price     = Price of the derivative

"""	
function distribution(mcProcess::BaseProcess,rfCurve::AbstractZeroRateCurve,mcConfig::MonteCarloConfiguration,T::num_) where { num_ <: Number }

	set_seed(mcConfig)
	S=simulate(mcProcess,rfCurve,mcConfig,T)

	return S[:,end];
end
