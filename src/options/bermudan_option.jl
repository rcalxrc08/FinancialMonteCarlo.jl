"""
Struct for Standard Bermudan Option

		bmOption=BermudanOption(T::AbstractArray{num1},K::num2,isCall::Bool=true) where {num1 <: Number,num2 <: Number}
	
Where:\n
		T	=	Time to maturity of the Option.
		K	=	Strike Price of the Option.
		isCall  = true for CALL, false for PUT.
"""
mutable struct BermudanOption{num1 <: Number,num2 <: Number}<:BermudanPayoff
	T::num1
	T_ex::Array{num1}
	K::num2
	isCall::Bool
	function BermudanOption(T::AbstractArray{num1},K::num2,isCall::Bool=true) where {num1 <: Number, num2 <: Number}
        if minimum(T) <= 0.0
            error("Time to Maturity must be positive")
        elseif K <= 0.0
            error("Strike Price must be positive")
        else
			@assert issorted(T);
            return new{num1,num2}(maximum(T),T,K,isCall)
        end
    end
end

export BermudanOption;

function payout(Sti::numtype_,bmPayoff::BermudanOption) where {numtype_<:Number}
	iscall=bmPayoff.isCall ? 1 : -1
	return ((Sti-bmPayoff.K)*iscall>0.0) ? (Sti-bmPayoff.K)*iscall : zero(numtype_);
end