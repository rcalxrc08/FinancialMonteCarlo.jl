type EuropeanOption<:EuropeanPayoff end

struct EUOptionData<:AbstractEuropeanOptionData
	T::Float64
	K::Float64
end

export EuropeanOption,EUOptionData;

"""
Payoff computation from MonteCarlo paths

		Payoff=payoff(S,euOptionData,EuropeanOption,isCall=true)
	
Where:\n
		S           = Paths of the Underlying.
		euOptionData  = Datas of the Option.
		EuropeanOption = Type of the Option
		isCall = true for Call Options, false for Put Options.

		Payoff      = payoff of the option.
```
"""
function payoff(S::Matrix{num},euOptionData::EUOptionData,spotData::equitySpotData,Payoff::EuropeanOption,isCall::Bool=true,T1::Float64=euOptionData.T) where{num<:Number}
	r=spotData.r;
	T=euOptionData.T;
	iscall=isCall?1:-1
	(Nsim,NStep)=size(S)
	NStep-=1;
	index1=round(Int,T/T1 * NStep)+1;
	S1=view(S,:,1:index1)
	ST=S1[:,end];
	K=euOptionData.K;
	f(ST::num)::num=(iscall*(ST-K)>0.0)?iscall*(ST-K):0.0;
	payoff2=f.(ST);
	
	return payoff2*exp(-r*T);
end