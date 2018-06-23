
struct BarrierOptionUpInData<:OptionData
	T::Float64
	K::Float64
	barrier::Float64
	isCall::Bool
	function BarrierOptionUpInData(T::Float64,K::Float64,barrier::Float64,isCall::Bool=true)
        if T <= 0.0
            error("Time to Maturity must be positive")
        elseif K <= 0.0
            error("Strike Price must be positive")
        elseif barrier <= 0.0
            error("Barrier must be positive")
        else
            return new(T,K,barrier,isCall)
        end
    end
end

export BarrierOptionUpInData;

"""
Payoff computation from MonteCarlo paths

		Payoff=payoff(S,barrierOptionData,)
	
Where:\n
		S           = Paths of the Underlying.
		barrierOptionData  = Datas of the Option.
		isCall = true for Call Options, false for Put Options.

		Payoff      = payoff of the option.
```
"""
function payoff(S::Matrix{num},barrierOptionData::BarrierOptionUpInData,spotData::equitySpotData,T1::Float64=barrierOptionData.T) where{num<:Number}
	iscall=barrierOptionData.isCall?1:-1
	r=spotData.r;
	T=barrierOptionData.T;
	(Nsim,NStep)=size(S)
	NStep-=1;
	K=barrierOptionData.K;
	D=barrierOptionData.barrier;
	index1=round(Int,T/T1 * NStep)+1;
	S1=view(S,:,1:index1)
	@inbounds f(S::Array{num})::num=(iscall*(S[end]-K)>0.0)&&(maximum(S)>D)?iscall*(S[end]-K):0.0;		
	@inbounds payoff2=[f(S1[i,:]) for i in 1:Nsim];
	
	return payoff2*exp(-r*T);
end