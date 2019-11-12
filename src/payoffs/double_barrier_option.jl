"""
Struct for Double Barrier Option

		dbOption=DoubleBarrierOption(T::num1,K::num2,D::num3,U::num4,isCall::Bool=true) where {num1 <: Number,num2 <: Number,num3 <: Number,num4 <: Number}
	
Where:\n
		T	=	Time to maturity of the Option.
		K	=	Strike Price of the Option.
		D	=	Low Barrier of the Option.
		U	=	Up Barrier of the Option.
		isCall  = true for CALL, false for PUT.
"""
mutable struct DoubleBarrierOption{num1 <: Number,num2 <: Number,num3 <: Number,num4 <: Number}<:BarrierPayoff
	T::num1
	K::num2
	D::num3
	U::num4
	isCall::Bool
	function DoubleBarrierOption(T::num1,K::num2,D::num3,U::num4,isCall::Bool=true) where {num1 <: Number,num2 <: Number,num3 <: Number,num4 <: Number}
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


function payoff(S::AbstractMatrix{num},doubleBarrierPayoff::DoubleBarrierOption,spotData::ZeroRateCurve,T1::num2=doubleBarrierPayoff.T) where{num <: Number,num2 <: Number}
	r=spotData.r;
	T=doubleBarrierPayoff.T;
	iscall=doubleBarrierPayoff.isCall ? 1 : -1
	(Nsim,NStep)=size(S)
	NStep-=1;
	K=doubleBarrierPayoff.K;
	D=doubleBarrierPayoff.D;
	U=doubleBarrierPayoff.U;
	index1=round(Int,T/T1 * NStep)+1;
	@inbounds f(S::AbstractArray{num})::num=(iscall*(S[end]-K)>0.0)&&(minimum(S)>D)&&(maximum(S)<U) ? iscall*(S[end]-K) : 0.0;		
	@inbounds payoff2=[f(S[i,1:index1]) for i in 1:Nsim];
	
	return payoff2*exp(-r*T);
end