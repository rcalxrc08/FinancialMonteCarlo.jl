"""
Struct for Binary European Option

		binOption=BinaryEuropeanOption(T::num1,K::num2,isCall::Bool=true) where {num1 <: Number,num2 <: Number}
	
Where:\n
		T	=	Time to maturity of the Option.
		K	=	Strike Price of the Option.
		isCall  = true for CALL, false for PUT.
"""
struct BinaryEuropeanOption{num1 <: Number,num2 <: Number}<:EuropeanPayoff
	T::num1
	K::num2
	isCall::Bool
	function BinaryEuropeanOption(T::num1,K::num2,isCall::Bool=true) where {num1 <: Number,num2 <: Number}
        if T <= 0.0
            error("Time to Maturity must be positive")
        elseif K <= 0.0
            error("Strike Price must be positive")
        else
            return new{num1,num2}(T,K,isCall)
        end
    end
end

export BinaryEuropeanOption;

function payout(ST::numtype_,euPayoff::BinaryEuropeanOption) where {numtype_<:Number}
	iscall=euPayoff.isCall ? 1 : -1
	return ((ST-euPayoff.K)*iscall>0.0) ? one(numtype_) : zero(numtype_);
end