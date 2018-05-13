type BinaryEuropeanOption<:EuropeanPayoff end

struct BinaryEuropeanOptionData<:AbstractEuropeanOptionData
	T::Float64
	K::Float64
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
function payoff(S::Matrix{num},euOptionData::BinaryEuropeanOption,spotData::equitySpotData,Payoff::BinaryEuropeanOption,isCall::Bool=true) where{num<:Number}
	r=spotData.r;
	T=euOptionData.T;
	iscall=isCall?1:-1
	ST=S[1:end,end];
	K=euOptionData.K;
	f(ST::num)::num=iscall*(ST-K)>0.0;
	payoff2=f.(ST);
	
	return payoff2*exp(-r*T);
end