
struct DoubleBarrierOption{num1,num2,num3,num4<:Number}<:BarrierPayoff
	T::num1
	K::num2
	D::num3
	U::num4
	isCall::Bool
	function DoubleBarrierOption(T::num1,K::num2,D::num3,U::num4,isCall::Bool=true) where {num1,num2,num3,num4<:Number}
        if T <= 0.0
            error("Time to Maturity must be positive")
        elseif K <= 0.0
            error("Strike Price must be positive")
        elseif D <= 0.0
            error("Low Barrier must be positive")
        elseif U <= 0.0
            error("High Barrier must be positive")
        else
            return new{num1,num2,num3,num4}(T,K,D,U,isCall)
        end
    end
end

export DoubleBarrierOption;

"""
Payoff computation from MonteCarlo paths

		Payoff=payoff(S,doubleBarrierPayoff,)
	
Where:\n
		S           = Paths of the Underlying.
		doubleBarrierPayoff  = Datas of the Option.

		Payoff      = payoff of the option.
```
"""
function payoff(S::AbstractMatrix{num},doubleBarrierPayoff::DoubleBarrierOption,spotData::equitySpotData,T1::num2=doubleBarrierPayoff.T) where{num,num2<:Number}
	r=spotData.r;
	T=doubleBarrierPayoff.T;
	iscall=doubleBarrierPayoff.isCall ? 1 : -1
	(Nsim,NStep)=size(S)
	NStep-=1;
	K=doubleBarrierPayoff.K;
	D=doubleBarrierPayoff.D;
	U=doubleBarrierPayoff.U;
	index1=round(Int,T/T1 * NStep)+1;
	@inbounds f(S::Array{num})::num=(iscall*(S[end]-K)>0.0)&&(minimum(S)>D)&&(maximum(S)<U) ? iscall*(S[end]-K) : 0.0;		
	@inbounds payoff2=[f(S[i,1:index1]) for i in 1:Nsim];
	
	return payoff2*exp(-r*T);
end