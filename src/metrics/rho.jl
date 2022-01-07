"""
General Interface for Rho computation

		Rho=rho(mcProcess,rfCurve,mcConfig,abstractPayoff,dr)
	
Where:\n
		mcProcess          = Process to be simulated.
		rfCurve  = Zero Rate Data.
		mcConfig = Basic properties of MonteCarlo simulation
		abstractPayoff = Payoff(s) to be priced
		dr [optional, default to 1e-7] = increment

		Rho     = Rho of the derivative

"""
function rho(mcProcess::BaseProcess, rfCurve::AbstractZeroRateCurve, mcConfig::MonteCarloConfiguration, abstractPayoff::AbstractPayoff, dr::Real = 1e-7)
    Price = pricer(mcProcess, rfCurve, mcConfig, abstractPayoff)
    rfCurve_1 = ZeroRate(rfCurve.r + dr)
    PriceUp = pricer(mcProcess, rfCurve_1, mcConfig, abstractPayoff)
    rho = (PriceUp - Price) / dr

    return rho
end

function rho(mcProcess::BaseProcess, rfCurve::AbstractZeroRateCurve, mcConfig::MonteCarloConfiguration, abstractPayoffs::Array{abstractPayoff_}, dr::Real = 1e-7) where {abstractPayoff_ <: AbstractPayoff}
    Prices = pricer(mcProcess, rfCurve, mcConfig, abstractPayoffs)
    rfCurve_1 = ZeroRate(rfCurve.r + dr)
    PricesUp = pricer(mcProcess, rfCurve_1, mcConfig, abstractPayoffs)
    rho = @. (PricesUp - Prices) / dr

    return rho
end

export rho;
