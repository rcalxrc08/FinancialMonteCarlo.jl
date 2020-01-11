"""
Struct for European Option

		euOption=EuropeanOptionND(T::num1,K::num2,isCall::Bool=true) where {num1 <: Number,num2 <: Number}
	
Where:\n
		T	=	Time to maturity of the Option.
		K	=	Strike Price of the Option.
		isCall  = true for CALL, false for PUT.
"""
mutable struct EuropeanOptionND{num1 <: Number ,num2<:Number}<:EuropeanBasketPayoff
	T::num1
	K::num2
	isCall::Bool
	function EuropeanOptionND(T::num1,K::num2,isCall::Bool=true) where {num1 <: Number , num2 <: Number}
        if T <= 0.0
            error("Time to Maturity must be positive")
        elseif K <= 0.0
            error("Strike Price must be positive")
        else
            return new{num1,num2}(T,K,isCall)
        end
    end
end

export EuropeanOptionND;


function payoff(S::Array{abstractMatrix},euPayoff::EuropeanOptionND,rfCurve::abstractZeroRateCurve,T1::num2=maturity(euPayoff)) where { abstractMatrix <: AbstractMatrix{num}, num2 <: Number} where {abstractZeroRateCurve <: AbstractZeroRateCurve, num <: Number}
	r=rfCurve.r;
	T=euPayoff.T;
	iscall=euPayoff.isCall ? 1 : -1
	(Nsim,NStep)=size(S[1])
	NStep-=1;
	index1=round(Int,T/T1 * NStep)+1;
	K=euPayoff.K;
	ST_all=[ sum(x_i[j,index1] for x_i in S) for j in 1:Nsim]
	payoff2=max.(iscall*(ST_all.-K),0.0);
	
	return payoff2*exp(-integral(r,T));
end