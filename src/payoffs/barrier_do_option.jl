type BarrierOptionDownOut<:BarrierPayoff end

struct BarrierOptionData<:OptionData
	T::Float64
	K::Float64
	barrier::Float64
	function BarrierOptionData(T::Float64,K::Float64,barrier::Float64)
        if T <= 0.0
            error("Time to Maturity must be positive")
        elseif K <= 0.0
            error("Strike Price must be positive")
        elseif barrier <= 0.0
            error("Barrier must be positive")
        else
            return new(T,K,barrier)
        end
    end
end

export BarrierOptionDownOut,BarrierOptionData;

"""
Payoff computation from MonteCarlo paths

		Payoff=payoff(S,barrierOptionData,BarrierOptionDownOut,isCall=true)

Where:\n
		S           = Paths of the Underlying.
		barrierOptionData  = Datas of the Option.
		BarrierOptionDownOut = Type of the Option
		isCall = true for Call Options, false for Put Options.

		Payoff      = payoff of the option.
```
"""
function payoff(S::Matrix{num},barrierOptionData::BarrierOptionData,spotData::equitySpotData,Payoff::BarrierOptionDownOut,isCall::Bool=true,T1::Float64=barrierOptionData.T) where{num<:Number}
	iscall=isCall?1:-1
	r=spotData.r;
	T=barrierOptionData.T;
	(Nsim,NStep)=size(S)
	NStep-=1;
	K=barrierOptionData.K;
	barrier=barrierOptionData.barrier;
	index1=round(Int,T/T1 * NStep)+1;
	S1=view(S,:,1:index1)
	@inbounds f(S::Array{num})::num=(iscall*(S[end]-K)>0.0)&&(minimum(S)>barrier)?iscall*(S[end]-K):0.0;
	@inbounds payoff2=[f(S1[i,:]) for i in 1:Nsim];

	return payoff2*exp(-r*T);
end
