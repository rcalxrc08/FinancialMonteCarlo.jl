
struct DoubleBarrierOption<:OptionData
	T::Float64
	K::Float64
	D::Float64
	U::Float64
	isCall::Bool
	function DoubleBarrierOption(T::Float64,K::Float64,D::Float64,U::Float64,isCall::Bool=true)
        if T <= 0.0
            error("Time to Maturity must be positive")
        elseif K <= 0.0
            error("Strike Price must be positive")
        elseif D <= 0.0
            error("Low Barrier must be positive")
        elseif U <= 0.0
            error("High Barrier must be positive")
        else
            return new(T,K,D,U,isCall)
        end
    end
end

export DoubleBarrierOption;

"""
Payoff computation from MonteCarlo paths

		Payoff=payoff(S,doubleBarrierOptionData,)
	
Where:\n
		S           = Paths of the Underlying.
		doubleBarrierOptionData  = Datas of the Option.

		Payoff      = payoff of the option.
```
"""
function payoff(S::Matrix{num},doubleBarrierOptionData::DoubleBarrierOption,spotData::equitySpotData,T1::Float64=doubleBarrierOptionData.T) where{num<:Number}
	r=spotData.r;
	T=doubleBarrierOptionData.T;
	iscall=doubleBarrierOptionData.isCall?1:-1
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