"""
Struct for Asian Fixed Strike Option

		asOption=AsianFixedStrikeOption(T::num1,K::num2,isCall::Bool=true) where {num1 <: Number,num2 <: Number}
	
Where:\n
		T	=	Time to maturity of the Option.
		K	=	Strike Price of the Option.
		isCall  = true for CALL, false for PUT.
"""
mutable struct AsianFixedStrikeOption{T1,T2<:Number}<:AsianPayoff
	T::T1
	K::T2
	isCall::Bool
	function AsianFixedStrikeOption(T::T1,K::T2,isCall::Bool=true) where {T1, T2 <: Number}
        if T <= 0.0
            error("Time to Maturity must be positive")
        elseif K <= 0.0
            error("Strike Price must be positive")
        else
            return new{T1,T2}(T,K,isCall)
        end
    end

end

export AsianFixedStrikeOption;

function payout(S::abstractArray,opt_::AsianFixedStrikeOption) where {abstractArray<:AbstractArray{num_}} where {num_<:Number}
	iscall=opt_.isCall ? 1 : -1
	zero_typed=zero(eltype(S))*opt_.K;
	return max(iscall*(mean(S)-opt_.K),zero_typed)
end