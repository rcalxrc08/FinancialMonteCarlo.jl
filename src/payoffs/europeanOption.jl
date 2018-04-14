type EuropeanOption<:PayoffMC end

struct EUOptionData<:OptionData
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
function payoff(S::Matrix{num},euOptionData::EUOptionData,Payoff::EuropeanOption,isCall::Bool=true) where{num<:Number}
	iscall=isCall?1:-1
	ST=S[1:end,end];
	K=euOptionData.K;
	f(ST::num)::num=(iscall*(ST-K)>0.0)?iscall*(ST-K):0.0;
	payoff2=f.(ST);
	
	return payoff2;
end