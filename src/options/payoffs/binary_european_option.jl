"""
Struct for Binary European Option

		binOption=BinaryEuropeanOption(T::num1,K::num2,isCall::Bool=true) where {num1 <: Number,num2 <: Number}
	
Where:

		T	    = Time to maturity of the Option.
		K	    = Strike Price of the Option.
		isCall  = true for CALL, false for PUT.
"""
struct BinaryEuropeanOption{num1 <: Number, num2 <: Number, numtype <: Number} <: EuropeanPayoff{numtype}
    T::num1
    K::num2
    isCall::Bool
    function BinaryEuropeanOption(T::num1, K::num2, isCall::Bool = true) where {num1 <: Number, num2 <: Number}
        ChainRulesCore.@ignore_derivatives @assert T > 0 "Time to Maturity must be positive"
        ChainRulesCore.@ignore_derivatives @assert K > 0 "Strike Price must be positive"
        zero_typed = ChainRulesCore.@ignore_derivatives zero(num1) + zero(num2)
        return new{num1, num2, typeof(zero_typed)}(T, K, isCall)
    end
end

export BinaryEuropeanOption;

function payout_untyped(ST::numtype_, euPayoff::BinaryEuropeanOption) where {numtype_ <: Number}
    iscall = ChainRulesCore.@ignore_derivatives ifelse(euPayoff.isCall, Int8(1), Int8(-1))
    res = ifelse((ST - 1) * iscall > 0, one(numtype_), zero(numtype_))
    return res
end

function payout(ST::numtype_, euPayoff::BinaryEuropeanOption) where {numtype_ <: Number}
    K = euPayoff.K
    return payout_untyped(ST / K, euPayoff)
end
