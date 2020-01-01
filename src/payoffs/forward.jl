"""
Class for Dispatching Forward Payoff

		forward=Forward(T::num) where {num<:Number}
	
Where:\n
		T	=	Time to maturity of the Forward.
"""
mutable struct Forward{num<:Number}<:EuropeanPayoff
	T::num
	function Forward(T::num) where {num<:Number}
        if T <= 0.0
            error("Time to Maturity must be positive")
        else
            return new{num}(T)
        end
    end
end

export Forward;

function payout(ST::numtype_,euPayoff::Forward) where {numtype_<:Number}
	return ST;
end

function payoff(S::AbstractMatrix{num},optionData::Forward,rfCurve::ZeroRateCurve,T1::num2=optionData.T) where{num <: Number,num2 <: Number}
	r=rfCurve.r;
	T=optionData.T;
	(Nsim,NStep)=size(S)
	NStep-=1;
	index1=round(Int,T/T1 * NStep)+1;
	@views ST=S[:,index1];
	
	return ST*exp(-integral(r,T));
end
