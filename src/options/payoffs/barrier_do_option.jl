"""
Struct for Barrier Down and Out Option

		barOption=BarrierOptionDownOut(T::num1,K::num2,barrier::num3,isCall::Bool=true) where {num1 <: Number,num2 <: Number,num3 <: Number}
	
Where:

		T	    = Time to maturity of the Option.
		K	    = Strike Price of the Option.
		barrier	= Down Barrier of the Option.
		isCall  = true for CALL, false for PUT.
"""
struct BarrierOptionDownOut{num1 <: Number, num2 <: Number, num3 <: Number, numtype <: Number} <: BarrierPayoff{numtype}
    T::num1
    K::num2
    barrier::num3
    isCall::Bool
    function BarrierOptionDownOut(T::num1, K::num2, barrier::num3, isCall::Bool = true) where {num1 <: Number, num2 <: Number, num3 <: Number}
        @assert T > 0 "Time to Maturity must be positive"
        @assert K > 0 "Strike Price must be positive"
        @assert barrier > 0 "Barrier must be positive"
        zero_typed = zero(num1) + zero(num2) + zero(num3)
        return new{num1, num2, num3, typeof(zero_typed)}(T, K, barrier, isCall)
    end
end

export BarrierOptionDownOut;

function payout(S::abstractArray, barrierPayoff::BarrierOptionDownOut) where {abstractArray <: AbstractArray{num_}} where {num_ <: Number}
    iscall = barrierPayoff.isCall ? 1 : -1
    zero_typed = zero(eltype(S)) * barrierPayoff.K * barrierPayoff.barrier
    return @views max(iscall * (S[end] - barrierPayoff.K), zero_typed) * (findfirst(x -> x < barrierPayoff.barrier, S) === nothing)
end
