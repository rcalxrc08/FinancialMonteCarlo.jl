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
function payoff(S::Matrix{num},asianFixedStrikeOptionData::AsianFixedStrikeOptionData,spotData::equitySpotData,Payoff::AsianFixedStrikeOption,isCall::Bool=true,T1::Float64=asianFixedStrikeOptionData.T) where{num<:Number}
	iscall=isCall?1:-1
	r=spotData.r;
	T=asianFixedStrikeOptionData.T;
	K=asianFixedStrikeOptionData.K;
	(Nsim,NStep)=size(S)
	NStep-=1;
	index1=round(Int,T/T1 * NStep)+1;
	S1=view(S,:,1:index1)
	@inbounds f(S::Array{num})::num=(iscall*(mean(S)-K)>0.0)?(iscall*(mean(S)-K)):0.0;		
	@inbounds payoff2=[f(S1[i,1:end]) for i in 1:Nsim];
	
	return payoff2*exp(-r*T);
end