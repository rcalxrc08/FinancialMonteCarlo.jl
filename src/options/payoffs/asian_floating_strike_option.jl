"""
Struct for Asian Floating Strike Option

		asOption=AsianFloatingStrikeOption(T::num1,isCall::Bool=true) where {num1 <: Number}
	
Where:

		T	    = Time to maturity of the Option.
		isCall  = true for CALL, false for PUT.
"""
struct AsianFloatingStrikeOption{num <: Number} <: AsianPayoff{num}
    T::num
    isCall::Bool
    function AsianFloatingStrikeOption(T::num, isCall::Bool = true) where {num <: Number}
        ChainRulesCore.@ignore_derivatives @assert T > 0 "Time to Maturity must be positive"
        return new{num}(T, isCall)
    end
end

export AsianFloatingStrikeOption;

function payout(S::abstractArray, opt_::AsianFloatingStrikeOption) where {abstractArray <: AbstractArray{num_}} where {num_ <: Number}
    iscall = ChainRulesCore.@ignore_derivatives ifelse(opt_.isCall, Int8(1), Int8(-1))
    zero_typed = ChainRulesCore.@ignore_derivatives zero(eltype(S))
    return max(iscall * (S[end] - mean(S)), zero_typed)
end
