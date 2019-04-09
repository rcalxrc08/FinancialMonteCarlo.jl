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

"""
Payoff computation from MonteCarlo paths

		Payoff=payoff(S,Forward,spotData,T1)

Where:\n
		S           = Paths of the Underlying.
		Forward  = Datas of the Forward.
		spotData  = Datas of the Spot.
		T1=Final Time of Spot Simulation (default equals Time to Maturity of the option)

		Payoff      = payoff of the Forward.
```
"""
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
