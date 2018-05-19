type AsianFloatingStrikeOption<:PayoffMC end

struct AsianFloatingStrikeOptionData<:OptionData
	T::Float64
end

export AsianFloatingStrikeOption,AsianFloatingStrikeOptionData;

"""
Payoff computation from MonteCarlo paths

		Payoff=payoff(S,asianFloatingStrikeOptionData,AsianFloatingStrikeOption,isCall=true)
	
Where:\n
		S           = Paths of the Underlying.
		asianFloatingStrikeOptionData  = Datas of the Option.
		AsianFloatingStrikeOption = Type of the Option
		isCall = true for Call Options, false for Put Options.

		Payoff      = payoff of the option.
```
"""
function payoff(S::Matrix{num},asianFloatingStrikeOptionData::AsianFloatingStrikeOptionData,spotData::equitySpotData,Payoff::AsianFloatingStrikeOption,isCall::Bool=true) where{num<:Number}
	iscall=isCall?1:-1
	r=spotData.r;
	T=asianFloatingStrikeOptionData.T;
	NsimTmp=length(S[1:end,end]);
	NStep=length(S[1,1:end])-1
	Nsim=length(S[1:end,1])
	dt=T1/NStep;
	index1=round(Int,T/dt);
	index1=index1>NStep?Nstep+1:index+1;
	S1=view(S,1:Nsim,1:index1)
	@inbounds f(S::Array{num})::num=(iscall*(S[end]-mean(S))>0.0)?(iscall*(S[end]-mean(S))):0.0;		
	@inbounds payoff2=[f(S1[i,1:end]) for i in 1:NsimTmp];
	
	return payoff2*exp(-r*T);
end