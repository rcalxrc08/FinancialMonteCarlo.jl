"""
Struct for European Option

		euOption=EuropeanOption(T::num1,K::num2,isCall::Bool=true) where {num1 <: Number,num2 <: Number}
	
Where:

		T	    = Time to maturity of the Option.
		K	    = Strike Price of the Option.
		isCall  = true for CALL, false for PUT.
"""
struct EuropeanOption{num1 <: Number, num2 <: Number, numtype <: Number} <: EuropeanPayoff{numtype}
    T::num1
    K::num2
    isCall::Bool
    function EuropeanOption(T::num1, K::num2, isCall::Bool = true) where {num1 <: Number, num2 <: Number}
        ChainRulesCore.@ignore_derivatives @assert T > 0 "Time to Maturity must be positive"
        ChainRulesCore.@ignore_derivatives @assert K > 0 "Strike Price must be positive"
        zero_typed = ChainRulesCore.@ignore_derivatives zero(num1) + zero(num2)
        return new{num1, num2, typeof(zero_typed)}(T, K, isCall)
    end
end

export EuropeanOption;

function payout_untyped(ST::numtype_, euPayoff::EuropeanOption) where {numtype_ <: Number}
    iscall = ChainRulesCore.@ignore_derivatives ifelse(euPayoff.isCall, Int8(1), Int8(-1))
    zero_typed = ChainRulesCore.@ignore_derivatives zero(ST)
    return max(iscall * (ST - 1), zero_typed)
end

function payout(ST::numtype_, euPayoff::EuropeanOption) where {numtype_ <: Number}
    K = euPayoff.K
    result = K * payout_untyped(ST / K, euPayoff)
    return result
end
