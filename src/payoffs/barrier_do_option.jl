
struct BarrierOptionDownOut<:BarrierPayoff
	T::Float64
	K::Float64
	barrier::Float64
	isCall::Bool
	function BarrierOptionDownOut(T::Float64,K::Float64,barrier::Float64,isCall::Bool=true)
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

export BarrierOptionDownOut;

"""
Payoff computation from MonteCarlo paths

		Payoff=payoff(S,barrierPayoff,)

Where:\n
		S           = Paths of the Underlying.
		barrierPayoff  = Datas of the Option.

		Payoff      = payoff of the option.
```
"""
function payoff(S::AbstractMatrix{num},barrierPayoff::BarrierOptionDownOut,spotData::equitySpotData,T1::Float64=barrierPayoff.T) where{num<:Number}
	iscall=barrierPayoff.isCall?1:-1
	r=spotData.r;
	T=barrierPayoff.T;
	(Nsim,NStep)=size(S)
	NStep-=1;
	K=barrierPayoff.K;
	barrier=barrierPayoff.barrier;
	index1=round(Int,T/T1 * NStep)+1;
	S1=view(S,:,1:index1)
	@inbounds f(S::Array{num})::num=(iscall*(S[end]-K)>0.0)&&(minimum(S)>barrier)?iscall*(S[end]-K):0.0;
	@inbounds payoff2=[f(S1[i,:]) for i in 1:Nsim];

	return payoff2*exp(-r*T);
end
