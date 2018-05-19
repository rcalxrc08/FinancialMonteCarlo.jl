type BarrierOptionDownOut<:PayoffMC end

struct BarrierOptionData<:OptionData
	T::Float64
	K::Float64
	D::Float64
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
	NsimTmp=length(S[1:end,end]);
	K=barrierOptionData.K;
	D=barrierOptionData.D;
	NStep=length(S[1,1:end])-1
	Nsim=length(S[1:end,1])
	dt=T1/NStep;
	index1=round(Int,T/dt);
	index1=index1>NStep?Nstep+1:index1+1;
	S1=view(S,1:Nsim,1:index1)
	@inbounds f(S::Array{num})::num=(iscall*(S[end]-K)>0.0)&&(minimum(S)>D)?iscall*(S[end]-K):0.0;		
	@inbounds payoff2=[f(S1[i,1:end]) for i in 1:NsimTmp];
	
	return payoff2*exp(-r*T);
end