"""
General Interface for distribution

		S=distribution(mcProcess,rfCurve,mcConfig,payoff_)
	
Where:\n
		mcProcess = Process to be simulated.
		rfCurve   = Zero Rate Data.
		mcConfig  = Basic properties of MonteCarlo simulation
		T         = Time of simulation

		S         = distribution

"""
function distribution(mcProcess::BaseProcess, rfCurve::AbstractZeroRateCurve, mcConfig::MonteCarloConfiguration, T::num_) where {num_ <: Number}
    set_seed!(mcConfig)
    S = simulate(mcProcess, rfCurve, mcConfig, T)

    return S[:, end]
end
