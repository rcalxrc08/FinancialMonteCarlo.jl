type BarrierOption<:PayoffMC end

struct BarrierOptionData<:OptionData
	T::Float64
	K::Float64
	D::Float64
end

export BarrierOption,BarrierOptionData;

"""
Payoff computation from MonteCarlo paths

		Payoff=payoff(S,barrierOptionData,BarrierOption,isCall=true)
	
Where:\n
		S           = Paths of the Underlying.
		barrierOptionData  = Datas of the Option.
		BarrierOption = Type of the Option
		isCall = true for Call Options, false for Put Options.

		Payoff      = payoff of the option.
```
"""
function payoff(S::Matrix{num},barrierOptionData::BarrierOptionData,spotData::equitySpotData,Payoff::BarrierOption,isCall::Bool=true) where{num<:Number}
	iscall=isCall?1:-1
	r=spotData.r;
	T=barrierOptionData.T;
	NsimTmp=length(S[1:end,end]);
	K=barrierOptionData.K;
	D=barrierOptionData.D;
	@inbounds f(S::Array{num})::num=(iscall*(S[end]-K)>0.0)&&(minimum(S)>D)?iscall*(S[end]-K):0.0;		
	@inbounds payoff2=[f(S[i,1:end]) for i in 1:NsimTmp];
	
	return payoff2*exp(-r*T);
end