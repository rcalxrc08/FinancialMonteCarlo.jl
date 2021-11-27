"""
Struct for Barrier Up and In Option

		barOption=BarrierOptionUpIn(T::num1,K::num2,barrier::num3,isCall::Bool=true) where {num1 <: Number,num2 <: Number,num3 <: Number}
	
Where:\n
		T	=	Time to maturity of the Option.
		K	=	Strike Price of the Option.
		barrier	=	Up Barrier of the Option.
		isCall  = true for CALL, false for PUT.
"""
struct BarrierOptionUpIn{num1 <: Number, num2 <: Number, num3 <: Number, numtype <: Number} <: BarrierPayoff{numtype}
    T::num1
    K::num2
    barrier::num3
    isCall::Bool
    function BarrierOptionUpIn(T::num1, K::num2, barrier::num3, isCall::Bool = true) where {num1 <: Number, num2 <: Number, num3 <: Number}
        if T <= 0.0
            error("Time to Maturity must be positive")
        elseif K <= 0.0
            error("Strike Price must be positive")
        elseif barrier <= 0.0
            error("Barrier must be positive")
        else
            zero_typed = zero(num1) + zero(num2) + zero(num3)
            return new{num1, num2, num3, typeof(zero_typed)}(T, K, barrier, isCall)
        end
    end
end

export BarrierOptionUpIn;

function payout(S::abstractArray, barrierPayoff::BarrierOptionUpIn) where {abstractArray <: AbstractArray{num_}} where {num_ <: Number}
    iscall = barrierPayoff.isCall ? 1 : -1
    zero_typed = zero(eltype(S)) * barrierPayoff.K * barrierPayoff.barrier
    return max(iscall * (S[end] - barrierPayoff.K), zero_typed) * (findfirst(x -> x > barrierPayoff.barrier, S) != nothing)
end
