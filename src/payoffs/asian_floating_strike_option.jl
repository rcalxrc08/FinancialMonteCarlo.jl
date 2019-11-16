"""
Struct for Asian Floating Strike Option

		asOption=AsianFloatingStrikeOption(T::num1,isCall::Bool=true) where {num1<:Number}
	
Where:\n
		T	=	Time to maturity of the Option.
		isCall  = true for CALL, false for PUT.
"""
mutable struct AsianFloatingStrikeOption{num<:Number}<:AsianPayoff
	T::num
	isCall::Bool
	function AsianFloatingStrikeOption(T::num,isCall::Bool=true) where {num<:Number}
        if T <= 0.0
            error("Time to Maturity must be positive")
        else
            return new{num}(T,isCall)
        end
    end
end

export AsianFloatingStrikeOption;


function payout(S::abstractArray,opt_::AsianFloatingStrikeOption) where {abstractArray<:AbstractArray{num_}} where {num_<:Number}
	iscall=opt_.isCall ? 1 : -1
	zero_typed=zero(eltype(S));
	return max(iscall*(S[end]-mean(S)),zero_typed)
end