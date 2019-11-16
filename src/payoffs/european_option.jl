"""
Struct for European Option

		euOption=EuropeanOption(T::num1,K::num2,isCall::Bool=true) where {num1 <: Number,num2 <: Number}
	
Where:\n
		T	=	Time to maturity of the Option.
		K	=	Strike Price of the Option.
		isCall  = true for CALL, false for PUT.
"""
struct EuropeanOption{num1 <: Number ,num2<:Number}<:EuropeanPayoff
	T::num1
	K::num2
	isCall::Bool
	function EuropeanOption(T::num1,K::num2,isCall::Bool=true) where {num1 <: Number , num2 <: Number}
        if T <= 0.0
            error("Time to Maturity must be positive")
        elseif K <= 0.0
            error("Strike Price must be positive")
        else
            return new{num1,num2}(T,K,isCall)
        end
    end
end

export EuropeanOption;


function payout(ST::numtype_,euPayoff::EuropeanOption) where {numtype_<:Number}
	iscall=euPayoff.isCall ? 1 : -1
	zero_typed=zero(ST)*euPayoff.K;
	return max.(iscall*(ST.-euPayoff.K),zero_typed);
end