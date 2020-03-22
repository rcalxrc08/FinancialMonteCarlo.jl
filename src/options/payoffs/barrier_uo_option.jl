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

function payout(S::abstractArray,barrierPayoff::BarrierOptionUpOut) where {abstractArray<:AbstractArray{num_}} where {num_<:Number}
	iscall=barrierPayoff.isCall ? 1 : -1
	zero_typed=zero(eltype(S))*barrierPayoff.K*barrierPayoff.barrier;
	return max(iscall*(S[end]-barrierPayoff.K),zero_typed)*(findfirst(x->x>barrierPayoff.barrier,S)==nothing);
end
