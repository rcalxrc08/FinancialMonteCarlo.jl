"""
Struct for Barrier Up and Out Option

		barOption=BarrierOptionUpOut(T::num1,K::num2,barrier::num3,isCall::Bool=true) where {num1 <: Number,num2 <: Number,num3 <: Number}
	
Where:\n
		T	=	Time to maturity of the Option.
		K	=	Strike Price of the Option.
		barrier	=	Up Barrier of the Option.
		isCall  = true for CALL, false for PUT.
"""
struct BarrierOptionUpOut{num1 <: Number,num2 <: Number,num3 <: Number}<:BarrierPayoff
	T::num1
	K::num2
	barrier::num3
	isCall::Bool
	function BarrierOptionUpOut(T::num1,K::num2,barrier::num3,isCall::Bool=true) where {num1 <: Number,num2 <: Number,num3 <: Number}
        if T <= 0.0
            error("Time to Maturity must be positive")
        elseif K <= 0.0
            error("Strike Price must be positive")
        elseif barrier <= 0.0
            error("Barrier must be positive")
        else
            return new{num1,num2,num3}(T,K,barrier,isCall)
        end
    end
end

export BarrierOptionUpOut;


function payoff(S::AbstractMatrix{num},barrierPayoff::BarrierOptionUpOut,spotData::ZeroRateCurve,T1::num2=barrierPayoff.T) where{num <: Number,num2 <: Number}
	iscall=barrierPayoff.isCall ? 1 : -1
	r=spotData.r;
	T=barrierPayoff.T;
	(Nsim,NStep)=size(S)
	NStep-=1;
	K=barrierPayoff.K;
	U=barrierPayoff.barrier;
	index1=round(Int,T/T1 * NStep)+1;
	@inbounds f(S::AbstractArray{num})::num=(iscall*(S[end]-K)>0.0)&&(maximum(S)<U) ? iscall*(S[end]-K) : 0.0;		
	@inbounds payoff2=[f(S[i,1:index1]) for i in 1:Nsim];
	
	return payoff2*exp(-r*T);
end