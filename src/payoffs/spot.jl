"""
Class for Dispatching Spot Payoff

		spot=Spot()
"""
mutable struct Spot <: EuropeanPayoff
	function Spot()
            return new()
    end
end

export Spot;


function payoff(S::AbstractMatrix{num},optionData::Spot,rfCurve::ZeroRateCurve,T1::num2=maturity(optionData)) where{num <: Number, num2 <: Number}

	S0_vec=S[:,1];
	
	return S0_vec;
end


function pricer(mcProcess::BaseProcess,rfCurve::AbstractZeroRateCurve,mcConfig::MonteCarloConfiguration,abstractPayoff::Spot)

	return mcProcess.underlying.S0;
end