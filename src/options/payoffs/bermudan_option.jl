"""
Struct for Standard Bermudan Option

		bmOption=BermudanOption(T::AbstractArray{num1},K::num2,isCall::Bool=true) where {num1 <: Number,num2 <: Number}
	
Where:\n
		T	=	Time to maturity of the Option.
		K	=	Strike Price of the Option.
		isCall  = true for CALL, false for PUT.
"""
mutable struct BermudanOption{num1 <: Number, num2 <: Number, numtype <: Number} <: BermudanPayoff{numtype}
    T::num1
    T_ex::Array{num1}
    K::num2
    isCall::Bool
    function BermudanOption(T::AbstractArray{num1}, K::num2, isCall::Bool = true) where {num1 <: Number, num2 <: Number}
        if minimum(T) <= 0.0
            error("Time to Maturity must be positive")
        elseif K <= 0.0
            error("Strike Price must be positive")
        else
            @assert issorted(T)
            zero_typed = zero(num1) + zero(num2)
            return new{num1, num2, typeof(zero_typed)}(T[end], T, K, isCall)
        end
    end
end

export BermudanOption;

function payout(Sti::numtype_, bmPayoff::BermudanOption) where {numtype_ <: Number}
    iscall = bmPayoff.isCall ? 1 : -1
    return ((Sti - bmPayoff.K) * iscall > 0.0) ? (Sti - bmPayoff.K) * iscall : zero(numtype_)
end
