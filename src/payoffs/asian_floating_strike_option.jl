type AsianFloatingStrikeOption<:AsianPayoff end

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
function payoff(S::Matrix{num},asianFloatingStrikeOptionData::AsianFloatingStrikeOptionData,spotData::equitySpotData,Payoff::AsianFloatingStrikeOption,isCall::Bool=true,T1::Float64=asianFloatingStrikeOptionData.T) where{num<:Number}
	iscall=isCall?1:-1
	r=spotData.r;
	T=asianFloatingStrikeOptionData.T;
	(Nsim,NStep)=size(S)
	NStep-=1;
	index1=round(Int,T/T1 * NStep)+1;
	S1=view(S,:,1:index1)
	@inbounds f(S::Array{num})::num=(iscall*(S[end]-mean(S))>0.0)?(iscall*(S[end]-mean(S))):0.0;		
	@inbounds payoff2=[f(S1[i,:]) for i in 1:Nsim];
	
	return payoff2*exp(-r*T);
end