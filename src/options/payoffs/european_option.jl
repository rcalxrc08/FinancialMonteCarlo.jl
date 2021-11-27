"""
Struct for European Option

		euOption=EuropeanOption(T::num1,K::num2,isCall::Bool=true) where {num1 <: Number,num2 <: Number}
	
Where:\n
		T	=	Time to maturity of the Option.
		K	=	Strike Price of the Option.
		isCall  = true for CALL, false for PUT.
"""
struct EuropeanOption{num1 <: Number, num2 <: Number, numtype <: Number} <: EuropeanPayoff{numtype}
    T::num1
    K::num2
    isCall::Bool
    function EuropeanOption(T::num1, K::num2, isCall::Bool = true) where {num1 <: Number, num2 <: Number}
        if T <= 0.0
            error("Time to Maturity must be positive")
        elseif K <= 0.0
            error("Strike Price must be positive")
        else
            zero_typed = zero(num1) + zero(num2)
            return new{num1, num2, typeof(zero_typed)}(T, K, isCall)
        end
    end
end

export EuropeanOption;

function payout(ST::numtype_, euPayoff::EuropeanOption) where {numtype_ <: Number}
    iscall = euPayoff.isCall ? 1 : -1
    zero_typed = zero(ST) * euPayoff.K
    return max(iscall * (ST - euPayoff.K), zero_typed)
end
