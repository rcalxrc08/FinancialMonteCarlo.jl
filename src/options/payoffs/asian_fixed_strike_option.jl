"""
Struct for Asian Fixed Strike Option

		asOption=AsianFixedStrikeOption(T::num1,K::num2,isCall::Bool=true) where {num1 <: Number,num2 <: Number}
	
Where:\n
		T	=	Time to maturity of the Option.
		K	=	Strike Price of the Option.
		isCall  = true for CALL, false for PUT.
"""
struct AsianFixedStrikeOption{num1 <: Number ,num2 <: Number, numtype <: Number} <: AsianPayoff{numtype}
	T::num1
	K::num2
	isCall::Bool
	function AsianFixedStrikeOption(T::num1,K::num2,isCall::Bool=true) where {num1 <: Number, num2 <: Number}
        if T <= 0.0
            error("Time to Maturity must be positive")
        elseif K <= 0.0
            error("Strike Price must be positive")
        else
			zero_typed=zero(num1)+zero(num2);
            return new{num1,num2,typeof(zero_typed)}(T,K,isCall)
        end
    end

end

export AsianFixedStrikeOption;

function payout(S::abstractArray,opt_::AsianFixedStrikeOption) where {abstractArray <: AbstractArray{num_}} where {num_ <: Number}
	iscall=opt_.isCall ? 1 : -1
	zero_typed=zero(eltype(S))*opt_.K;
	return max(iscall*(mean(S)-opt_.K),zero_typed)
end