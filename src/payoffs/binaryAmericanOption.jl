type BinaryAmericanOption<:PayoffMC end

struct AMOptionData<:AbstractEuropeanOptionData
	T::Float64
	K::Float64
end

export BinaryAmericanOption,AMOptionData;

"""
Payoff computation from MonteCarlo paths

		Payoff=payoff(S,amOptionData,BinaryAmericanOption,isCall=true)
	
Where:\n
		S           = Paths of the Underlying.
		amOptionData  = Datas of the Option.
		BinaryAmericanOption = Type of the Option
		isCall = true for Call Options, false for price Options.

		Payoff      = payoff of the option.
```
"""
function payoff(S::Matrix{num},amOptionData::AMOptionData,spotData::equitySpotData,Payoff::BinaryAmericanOption,isCall::Bool=true,T1::Float64=amOptionData.T) where{num<:Number}
	iscall=isCall?1:-1
	Nsim=length(S[1:end,1]);
	Nstep=length(S[1,1:end])
	K=amOptionData.K;
	T=amOptionData.T;
	NStep=length(S[1,1:end])-1
	Nsim=length(S[1:end,1])
	dt=T1/NStep;
	index1=round(Int,T/dt);
	index1=index1>NStep?Nstep+1:index1+1;
	S1=view(S,1:Nsim,1:index1)
	phi(Sti::Number)::Number=((Sti-K)*iscall>0.0)?1.0:0.0;
	
	payoff1=payoff(collect(S1),spotData,GeneralAmericanOption(),phi,T);
	
	return payoff1;
end
