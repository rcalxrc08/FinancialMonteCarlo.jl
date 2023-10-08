"""
Struct for Standard American Option

		amOption=AmericanOption(T::num1,K::num2,isCall::Bool=true) where {num1 <: Number,num2 <: Number}
	
Where:

		T	    = Time to maturity of the Option.
		K	    = Strike Price of the Option.
		isCall  = true for CALL, false for PUT.
"""
struct AmericanOption{num1 <: Number, num2 <: Number, numtype <: Number} <: AmericanPayoff{numtype}
    T::num1
    K::num2
    isCall::Bool
    function AmericanOption(T::num1, K::num2, isCall::Bool = true) where {num1 <: Number, num2 <: Number}
        ChainRulesCore.@ignore_derivatives @assert T > 0 "Time to Maturity must be positive"
        ChainRulesCore.@ignore_derivatives @assert K > 0 "Strike Price must be positive"
        zero_typed = ChainRulesCore.@ignore_derivatives zero(num1) + zero(num2)
        return new{num1, num2, typeof(zero_typed)}(T, K, isCall)
    end
end

export AmericanOption;

function payout(Sti::numtype_, amPayoff::AmericanOption) where {numtype_ <: Number}
    iscall = ChainRulesCore.@ignore_derivatives ifelse(amPayoff.isCall, Int8(1), Int8(-1))
    zero_typed = ChainRulesCore.@ignore_derivatives zero(Sti) * amPayoff.K
    return max(iscall * (Sti - amPayoff.K), zero_typed)
end
