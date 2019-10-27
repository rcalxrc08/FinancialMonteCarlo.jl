"""
Class for Dispatching Spot Payoff

		forward=Spot(T::num) where {num<:Number}
	
Where:\n
		T	=	Time to maturity of the Spot.
"""
struct Spot <: EuropeanPayoff
	function Spot()
            return new()
    end
end

export Spot;


function payoff(S::AbstractMatrix{num},optionData::Spot,spotData::equitySpotData,T1::num2=optionData.T) where{num <: Number,num2 <: Number}

	S0_vec=S[:,1];
	
	return S0_vec;
end
