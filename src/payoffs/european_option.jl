"""
Struct for European Option

		euOption=EuropeanOption(T::num1,K::num2,isCall::Bool=true) where {num1 <: Number,num2 <: Number}
	
Where:\n
		T	=	Time to maturity of the Option.
		K	=	Strike Price of the Option.
		isCall  = true for CALL, false for PUT.
"""
mutable struct EuropeanOption{num1 <: Number ,num2<:Number}<:EuropeanPayoff
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


function payoff(S::AbstractMatrix{num},euPayoff::EuropeanOption,spotData::ZeroRateCurve,T1::num2=euPayoff.T) where{ num <: Number, num2 <: Number}
	r=spotData.r;
	T=euPayoff.T;
	iscall=euPayoff.isCall ? 1 : -1
	(Nsim,NStep)=size(S)
	NStep-=1;
	index1=round(UInt,T/T1 * NStep)+1;
	@views ST=S[:,index1];
	K=euPayoff.K;
	zero_typed=zero(eltype(ST))*K;
	payoff2=max.(iscall*(ST.-K),zero_typed);
	
	return payoff2*exp(-r*T);
end