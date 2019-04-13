"""
Class for Dispatching Forward Payoff

		forward=Forward(T::num) where {num<:Number}
	
Where:\n
		T	=	Time to maturity of the Forward.
"""
struct BinaryEuropeanOption{num1,num2<:Number}<:EuropeanPayoff
	T::num1
	K::num2
	isCall::Bool
	function BinaryEuropeanOption(T::num1,K::num2,isCall::Bool=true) where {num1,num2<:Number}
        if T <= 0.0
            error("Time to Maturity must be positive")
        elseif K <= 0.0
            error("Strike Price must be positive")
        else
            return new{num1,num2}(T,K,isCall)
        end
    end
end

export BinaryEuropeanOption;

"""
Payoff computation from MonteCarlo paths

		Payoff=payoff(S,EuropeanOption,spotData,T1)

Where:\n
		S           = Paths of the Underlying.
		EuropeanOption  = Datas of the Option.
		spotData  = Datas of the Spot.
		T1=Final Time of Spot Simulation (default equals Time to Maturity of the option)

		Payoff      = payoff of the Option.
```
"""
function payoff(S::AbstractMatrix{num},euPayoff::BinaryEuropeanOption,spotData::equitySpotData,T1::num2=euPayoff.T) where{num,num2<:Number}
	r=spotData.r;
	T=euPayoff.T;
	iscall=euPayoff.isCall ? 1 : -1
	K=euPayoff.K;
	(Nsim,NStep)=size(S)
	NStep-=1;
	index1=round(Int,T/T1 * NStep)+1;
	S1=view(S,:,1:index1)
	ST=S1[:,end];
	f(ST::num)::num=iscall*(ST-K)>0.0;
	payoff2=f.(ST);

	return payoff2*exp(-r*T);
end
