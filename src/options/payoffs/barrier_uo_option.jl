"""
Struct for Barrier Up and Out Option

		barOption=BarrierOptionUpOut(T::num1,K::num2,barrier::num3,isCall::Bool=true) where {num1 <: Number,num2 <: Number,num3 <: Number}
	
Where:

		T	    = Time to maturity of the Option.
		K	    = Strike Price of the Option.
		barrier	= Up Barrier of the Option.
		isCall  = true for CALL, false for PUT.
"""
struct BarrierOptionUpOut{num1 <: Number, num2 <: Number, num3 <: Number, numtype <: Number} <: BarrierPayoff{numtype}
    T::num1
    K::num2
    barrier::num3
    isCall::Bool
    function BarrierOptionUpOut(T::num1, K::num2, barrier::num3, isCall::Bool = true) where {num1 <: Number, num2 <: Number, num3 <: Number}
        ChainRulesCore.@ignore_derivatives @assert T > 0 "Time to Maturity must be positive"
        ChainRulesCore.@ignore_derivatives @assert K > 0 "Strike Price must be positive"
        ChainRulesCore.@ignore_derivatives @assert barrier > 0 "Barrier must be positive"
        zero_typed = ChainRulesCore.@ignore_derivatives zero(num1) + zero(num2) + zero(num3)
        return new{num1, num2, num3, typeof(zero_typed)}(T, K, barrier, isCall)
    end
end

export BarrierOptionUpOut;

function payout(S::abstractArray, barrierPayoff::BarrierOptionUpOut) where {abstractArray <: AbstractArray{num_}} where {num_ <: Number}
    iscall = ChainRulesCore.@ignore_derivatives ifelse(barrierPayoff.isCall, Int8(1), Int8(-1))
    zero_typed = ChainRulesCore.@ignore_derivatives zero(eltype(S)) * barrierPayoff.K * barrierPayoff.barrier
    return max(iscall * (S[end] - barrierPayoff.K), zero_typed) * (findfirst(x -> x > barrierPayoff.barrier, S) === nothing)
end
