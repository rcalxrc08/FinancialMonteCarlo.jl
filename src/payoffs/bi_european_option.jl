"""
Struct for European Option

		euOption=EuropeanOption2D(T::num1,K::num2,isCall::Bool=true) where {num1 <: Number,num2 <: Number}
	
Where:\n
		T	=	Time to maturity of the Option.
		K	=	Strike Price of the Option.
		isCall  = true for CALL, false for PUT.
"""
struct EuropeanOption2D{num1 <: Number ,num2<:Number}<:EuropeanPayoff
	T::num1
	K::num2
	isCall::Bool
	function EuropeanOption2D(T::num1,K::num2,isCall::Bool=true) where {num1 <: Number , num2 <: Number}
        if T <= 0.0
            error("Time to Maturity must be positive")
        elseif K <= 0.0
            error("Strike Price must be positive")
        else
            return new{num1,num2}(T,K,isCall)
        end
    end
end

export EuropeanOption2D;


function payoff(S::Tuple{abstractMatrix,abstractMatrix},euPayoff::EuropeanOption2D,spotData::ZeroRateCurve,T1::num2=euPayoff.T) where { abstractMatrix <: AbstractMatrix{num}, num2 <: Number} where { num <: Number}
	r=spotData.r;
	T=euPayoff.T;
	iscall=euPayoff.isCall ? 1 : -1
	(Nsim,NStep)=size(S[1])
	NStep-=1;
	index1=round(Int,T/T1 * NStep)+1;
	ST=S[1][:,index1]+S[2][:,index1];
	K=euPayoff.K;
	payoff2=max.(iscall*(ST.-K),0.0);
	
	return payoff2*exp(-r*T);
end