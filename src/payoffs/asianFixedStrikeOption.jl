type AsianFixedStrikeOption<:PayoffMC end

struct AsianFixedStrikeOptionData<:OptionData
	T::Float64
	K::Float64
end

export AsianFixedStrikeOption,AsianFixedStrikeOptionData;

"""
Payoff computation from MonteCarlo paths

		Payoff=payoff(S,asianFixedStrikeOptionData,AsianFixedStrikeOption,isCall=true)
	
Where:\n
		S           = Paths of the Underlying.
		asianFixedStrikeOptionData  = Datas of the Option.
		AsianFixedStrikeOption = Type of the Option
		isCall = true for Call Options, false for Put Options.

		Payoff      = payoff of the option.
```
"""
function payoff(S::Matrix{num},asianFixedStrikeOptionData::AsianFixedStrikeOptionData,spotData::equitySpotData,Payoff::AsianFixedStrikeOption,isCall::Bool=true) where{num<:Number}
	iscall=isCall?1:-1
	r=spotData.r;
	T=asianFixedStrikeOptionData.T;
	K=asianFixedStrikeOptionData.K;
	NStep=length(S[1,1:end])-1
	Nsim=length(S[1:end,1])
	dt=T1/NStep;
	index1=round(Int,T/dt);
	index1=index1>NStep?Nstep+1:index+1;
	S1=view(S,1:Nsim,1:index1)
	NsimTmp=length(S[1:end,end]);
	@inbounds f(S::Array{num})::num=(iscall*(mean(S)-K)>0.0)?(iscall*(mean(S)-K)):0.0;		
	@inbounds payoff2=[f(S1[i,1:end]) for i in 1:NsimTmp];
	
	return payoff2*exp(-r*T);
end