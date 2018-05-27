type BinaryEuropeanOption<:EuropeanPayoff end

struct BinaryEuropeanOptionData<:AbstractEuropeanOptionData
	T::Float64
	K::Float64
	BinaryEuropeanOptionData(T,K)=ifelse((T<=0.0),error("Cannot have negative Time to Maturity"),ifelse((K<=0.0),error("Cannot have negative Strike Price"),new(T,K)))
end

export BinaryEuropeanOption,BinaryEuropeanOptionData;

"""
Payoff computation from MonteCarlo paths

		Payoff=payoff(S,euOptionData,EuropeanOption,isCall=true)
	
Where:\n
		S           = Paths of the Underlying.
		euOptionData  = Datas of the Option.
		BinaryEuropeanOption = Type of the Option
		isCall = true for Call Options, false for Put Options.

		Payoff      = payoff of the option.
```
"""
function payoff(S::Matrix{num},euOptionData::BinaryEuropeanOption,spotData::equitySpotData,Payoff::BinaryEuropeanOption,isCall::Bool=true,T1::Float64=euOptionData.T) where{num<:Number}
	r=spotData.r;
	T=euOptionData.T;
	iscall=isCall?1:-1
	K=euOptionData.K;
	(Nsim,NStep)=size(S)
	NStep-=1;
	index1=round(Int,T/T1 * NStep)+1;
	S1=view(S,:,1:index1)
	ST=S1[:,end];
	f(ST::num)::num=iscall*(ST-K)>0.0;
	payoff2=f.(ST);
	
	return payoff2*exp(-r*T);
end