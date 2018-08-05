
struct AsianFloatingStrikeOption{num<:Number}<:AsianPayoff
	T::num
	isCall::Bool
	function AsianFloatingStrikeOption(T::num,isCall::Bool=true) where {num<:Number}
        if T <= 0.0
            error("Time to Maturity must be positive")
        else
            return new{num}(T,isCall)
        end
    end
end

export AsianFloatingStrikeOption;

"""
Payoff computation from MonteCarlo paths

		Payoff=payoff(S,asianFloatingStrikePayoff,)

Where:\n
		S           = Paths of the Underlying.
		asianFloatingStrikePayoff  = Datas of the Option.

		Payoff      = payoff of the option.
```
"""
function payoff(S::AbstractMatrix{num},asianFloatingStrikePayoff::AsianFloatingStrikeOption,spotData::equitySpotData,T1::num2=asianFloatingStrikePayoff.T) where{num,num2<:Number}
	iscall=asianFloatingStrikePayoff.isCall ? 1 : -1
	r=spotData.r;
	T=asianFloatingStrikePayoff.T;
	(Nsim,NStep)=size(S)
	NStep-=1;
	index1=round(Int,T/T1 * NStep)+1;
	S1=view(S,:,1:index1)
	@inbounds f(S::Array{num})::num=(iscall*(S[end]-mean(S))>0.0) ? (iscall*(S[end]-mean(S))) : 0.0;
	@inbounds payoff2=[f(S1[i,:]) for i in 1:Nsim];

	return payoff2*exp(-r*T);
end
