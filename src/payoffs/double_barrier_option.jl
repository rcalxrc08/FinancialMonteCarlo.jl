type DoubleBarrierOption<:BarrierPayoff end

struct DoubleBarrierOptionData<:OptionData
	T::Float64
	K::Float64
	D::Float64
	U::Float64
	DoubleBarrierOptionData(T,K,D,U)=ifelse((T<=0.0),error("Cannot have negative Time to Maturity"),ifelse((K<=0.0),error("Cannot have negative Strike Price"),ifelse((D<=0.0),error("Cannot have negative Low Barrier Price"),ifelse((U<=0.0),error("Cannot have negative High Barrier Price"),new(T,K,D,U)))))
end

export DoubleBarrierOption,DoubleBarrierOptionData;

"""
Payoff computation from MonteCarlo paths

		Payoff=payoff(S,doubleBarrierOptionData,DoubleBarrierOption,isCall=true)
	
Where:\n
		S           = Paths of the Underlying.
		doubleBarrierOptionData  = Datas of the Option.
		DoubleBarrierOption = Type of the Option
		isCall = true for Call Options, false for Put Options.

		Payoff      = payoff of the option.
```
"""
function payoff(S::Matrix{num},doubleBarrierOptionData::DoubleBarrierOptionData,spotData::equitySpotData,Payoff::DoubleBarrierOption,isCall::Bool=true,T1::Float64=doubleBarrierOptionData.T) where{num<:Number}
	r=spotData.r;
	T=doubleBarrierOptionData.T;
	iscall=isCall?1:-1
	(Nsim,NStep)=size(S)
	NStep-=1;
	K=doubleBarrierOptionData.K;
	D=doubleBarrierOptionData.D;
	U=doubleBarrierOptionData.U;
	index1=round(Int,T/T1 * NStep)+1;
	@inbounds f(S::Array{num})::num=(iscall*(S[end]-K)>0.0)&&(minimum(S)>D)&&(maximum(S)<U)?iscall*(S[end]-K):0.0;		
	@inbounds payoff2=[f(S[i,1:index1]) for i in 1:Nsim];
	
	return payoff2*exp(-r*T);
end