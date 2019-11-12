"""
Struct for Binary European Option

		binOption=BinaryEuropeanOption(T::num1,K::num2,isCall::Bool=true) where {num1 <: Number,num2 <: Number}
	
Where:\n
		T	=	Time to maturity of the Option.
		K	=	Strike Price of the Option.
		isCall  = true for CALL, false for PUT.
"""
mutable struct BinaryEuropeanOption{num1 <: Number,num2 <: Number}<:EuropeanPayoff
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


function payoff(S::AbstractMatrix{num},euPayoff::BinaryEuropeanOption,spotData::ZeroRateCurve,T1::num2=euPayoff.T) where{num <: Number,num2 <: Number}
	r=spotData.r;
	T=euPayoff.T;
	iscall=euPayoff.isCall ? 1 : -1
	K=euPayoff.K;
	(Nsim,NStep)=size(S)
	NStep-=1;
	index1=round(Int,T/T1 * NStep)+1;
	ST=S[:,index1];
	f(ST::num)::num=iscall*(ST-K)>0.0;
	payoff2=f.(ST);

	return payoff2*exp(-r*T);
end
