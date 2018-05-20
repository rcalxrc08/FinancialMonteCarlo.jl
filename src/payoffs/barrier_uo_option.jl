type BarrierOptionUpOut<:BarrierPayoff end

export BarrierOptionUpOut;

"""
Payoff computation from MonteCarlo paths

		Payoff=payoff(S,barrierOptionData,BarrierOptionUpOut,isCall=true)
	
Where:\n
		S           = Paths of the Underlying.
		barrierOptionData  = Datas of the Option.
		BarrierOptionUpOut = Type of the Option
		isCall = true for Call Options, false for Put Options.

		Payoff      = payoff of the option.
```
"""
function payoff(S::Matrix{num},barrierOptionData::BarrierOptionData,spotData::equitySpotData,Payoff::BarrierOptionUpOut,isCall::Bool=true,T1::Float64=barrierOptionData.T) where{num<:Number}
	iscall=isCall?1:-1
	r=spotData.r;
	T=barrierOptionData.T;
	(Nsim,NStep)=size(S)
	NStep-=1;
	K=barrierOptionData.K;
	U=barrierOptionData.barrier;
	index1=round(Int,T/T1 * NStep)+1;
	S1=view(S,:,1:index1)
	@inbounds f(S::Array{num})::num=(iscall*(S[end]-K)>0.0)&&(maximum(S)<U)?iscall*(S[end]-K):0.0;		
	@inbounds payoff2=[f(S1[i,:]) for i in 1:Nsim];
	
	return payoff2*exp(-r*T);
end