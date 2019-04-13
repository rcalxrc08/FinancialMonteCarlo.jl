"""
Class for Dispatching Forward Payoff

		forward=Forward(T::num) where {num<:Number}
	
Where:\n
		T	=	Time to maturity of the Forward.
"""
struct Forward{num<:Number}<:EuropeanPayoff
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


function payoff(S::AbstractMatrix{num},optionData::Forward,spotData::equitySpotData,T1::num2=optionData.T) where{num,num2<:Number}
	r=spotData.r;
	T=optionData.T;
	(Nsim,NStep)=size(S)
	NStep-=1;
	index1=round(Int,T/T1 * NStep)+1;
	S1=view(S,:,1:index1)
	ST=S1[:,end];
	
	return ST*exp(-r*T);
end
