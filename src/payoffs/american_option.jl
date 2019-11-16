"""
Struct for Standard American Option

		amOption=AmericanOption(T::num1,K::num2,isCall::Bool=true) where {num1 <: Number,num2 <: Number}
	
Where:\n
		T	=	Time to maturity of the Option.
		K	=	Strike Price of the Option.
		isCall  = true for CALL, false for PUT.
"""
mutable struct AmericanOption{num1 <: Number,num2 <: Number}<:AmericanPayoff
	T::num1
	K::num2
	isCall::Bool
	function AmericanOption(T::num1,K::num2,isCall::Bool=true) where {num1 <: Number, num2 <: Number}
        if T <= 0.0
            error("Time to Maturity must be positive")
        elseif K <= 0.0
            error("Strike Price must be positive")
        else
            return new{num1,num2}(T,K,isCall)
        end
    end
end

export AmericanOption;

function payout(Sti::numtype_,amPayoff::AmericanOption) where {numtype_<:Number}
	iscall=amPayoff.isCall ? 1 : -1
	return ((Sti-amPayoff.K)*iscall>0.0) ? (Sti-amPayoff.K)*iscall : zero(numtype_);
end

function payoff(S::AbstractMatrix{numtype},amPayoff::AmericanOption,spotData::ZeroRateCurve,T1::num2=amPayoff.T) where{numtype <: Number,num2 <: Number}
	iscall=amPayoff.isCall ? 1 : -1
	K=amPayoff.K;
	T=amPayoff.T;
	isDualZero=typeof(zero(eltype(S))*K)
	
	#phi_(Sti::numtype_) where {numtype_<:Number}=((Sti-K)*iscall>0.0) ? (Sti-K)*iscall : zero(isDualZero);
	(Nsim,NStep)=size(S)
	NStep-=1;
	index1=round(Int,T/T1 * NStep)+1;#round(Int,T/T1 * NStep)
	payoff1=payoff(collect(S[:,1:index1]),GeneralAmericanOption(amPayoff),spotData,T);
	
	return payoff1;
end
