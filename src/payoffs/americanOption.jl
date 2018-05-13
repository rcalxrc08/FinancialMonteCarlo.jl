type AmericanOption<:PayoffMC end

struct AMOptionData<:AbstractEuropeanOptionData
	T::Float64
	K::Float64
end

export AmericanOption,AMOptionData;

"""
Payoff computation from MonteCarlo paths

		Payoff=payoff(S,amOptionData,AmericanOption,isCall=true)
	
Where:\n
		S           = Paths of the Underlying.
		amOptionData  = Datas of the Option.
		AmericanOption = Type of the Option
		isCall = true for Call Options, false for price Options.

		Payoff      = payoff of the option.
```
"""
function payoff(S::Matrix{num},amOptionData::AMOptionData,spotData::equitySpotData,Payoff::AmericanOption,isCall::Bool=true) where{num<:Number}
	iscall=isCall?1:-1
	Nsim=length(S[1:end,1]);
	Nstep=length(S[1,1:end])
	K=amOptionData.K;
	T=amOptionData.T;
	dt=T/Nstep
	phi(Sti::Number)::Number=((Sti-K)*iscall>0.0)?(Sti-K)*iscall:0.0;
	
	payoff1=payoff(S,spotData,GeneralAmericanOption(),phi,T);
	
	return payoff1;
end
