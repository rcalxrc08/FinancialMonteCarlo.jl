"""
Struct for Double Barrier Option

		dbOption=DoubleBarrierOption(T::num1,K::num2,D::num3,U::num4,isCall::Bool=true) where {num1 <: Number,num2 <: Number,num3 <: Number,num4 <: Number}
	
Where:

		T	    = Time to maturity of the Option.
		K	    = Strike Price of the Option.
		D	    = Low Barrier of the Option.
		U	    = Up Barrier of the Option.
		isCall  = true for CALL, false for PUT.
"""
struct DoubleBarrierOption{num1 <: Number, num2 <: Number, num3 <: Number, num4 <: Number, numtype <: Number} <: BarrierPayoff{numtype}
    T::num1
    K::num2
    D::num3
    U::num4
    isCall::Bool
    function DoubleBarrierOption(T::num1, K::num2, D::num3, U::num4, isCall::Bool = true) where {num1 <: Number, num2 <: Number, num3 <: Number, num4 <: Number}
        ChainRulesCore.@ignore_derivatives @assert T > 0 "Time to Maturity must be positive"
        ChainRulesCore.@ignore_derivatives @assert K > 0 "Strike Price must be positive"
        ChainRulesCore.@ignore_derivatives @assert D > 0 "Low Barrier must be positive"
        ChainRulesCore.@ignore_derivatives @assert U > D "High Barrier must be positive and greater then Low barrier"
        zero_typed = ChainRulesCore.@ignore_derivatives zero(num1) + zero(num2) + zero(num3) + zero(num4)
        return new{num1, num2, num3, num4, typeof(zero_typed)}(T, K, D, U, isCall)
    end
end

export DoubleBarrierOption;

function payout(S::abstractArray, barrierPayoff::DoubleBarrierOption) where {abstractArray <: AbstractArray{num_}} where {num_ <: Number}
    iscall = ChainRulesCore.@ignore_derivatives ifelse(barrierPayoff.isCall, Int8(1), Int8(-1))
    zero_typed = ChainRulesCore.@ignore_derivatives zero(eltype(S)) * barrierPayoff.K * barrierPayoff.U * barrierPayoff.D
    return max(iscall * (S[end] - barrierPayoff.K), zero_typed) * (findfirst(x -> x < barrierPayoff.D, S) === nothing) * (findfirst(x -> x > barrierPayoff.U, S) === nothing)
end
