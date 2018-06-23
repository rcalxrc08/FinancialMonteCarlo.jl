
struct AsianFloatingStrikeOptionData<:OptionData
	T::Float64
	isCall::Bool
	function AsianFloatingStrikeOptionData(T::Float64,isCall::Bool=true)
        if T <= 0.0
            error("Time to Maturity must be positive")
        else
            return new(T,isCall)
        end
    end
end

export AsianFloatingStrikeOptionData;

"""
Payoff computation from MonteCarlo paths

		Payoff=payoff(S,asianFloatingStrikeOptionData,)

Where:\n
		S           = Paths of the Underlying.
		asianFloatingStrikeOptionData  = Datas of the Option.

		Payoff      = payoff of the option.
```
"""
function payoff(S::Matrix{num},asianFloatingStrikeOptionData::AsianFloatingStrikeOptionData,spotData::equitySpotData,T1::Float64=asianFloatingStrikeOptionData.T) where{num<:Number}
	iscall=asianFloatingStrikeOptionData.isCall?1:-1
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
