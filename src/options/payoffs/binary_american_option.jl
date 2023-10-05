"""
Struct for Binary American Option

		binAmOption=BinaryAmericanOption(T::num1,K::num2,isCall::Bool=true) where {num1 <: Number,num2 <: Number}
	
Where:

		T	    = Time to maturity of the Option.
		K	    = Strike Price of the Option.
		isCall  = true for CALL, false for PUT.
"""
struct BinaryAmericanOption{num1 <: Number, num2 <: Number, numtype <: Number} <: AmericanPayoff{numtype}
    T::num1
    K::num2
    isCall::Bool
    function BinaryAmericanOption(T::num1, K::num2, isCall::Bool = true) where {num1 <: Number, num2 <: Number}
        @assert T > 0 "Time to Maturity must be positive"
        @assert K > 0 "Strike Price must be positive"
        zero_typed = zero(num1) + zero(num2)
        return new{num1, num2, typeof(zero_typed)}(T, K, isCall)
    end
end

export BinaryAmericanOption;

function payout(Sti::numtype_, amPayoff::BinaryAmericanOption) where {numtype_ <: Number}
    iscall = amPayoff.isCall ? 1 : -1
    return ((Sti - amPayoff.K) * iscall > 0.0) ? one(numtype_) : zero(numtype_)
end
