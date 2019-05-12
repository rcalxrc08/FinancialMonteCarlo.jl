"""
Struct for European Option

		euOption=EuropeanOption(T::num1,K::num2,isCall::Bool=true) where {num1,num2<:Number}
	
Where:\n
		T	=	Time to maturity of the Option.
		K	=	Strike Price of the Option.
		isCall  = true for CALL, false for PUT.
"""
struct EuropeanOption{num1,num2<:Number}<:EuropeanPayoff
	T::num1
	K::num2
	isCall::Bool
	function EuropeanOption(T::num1,K::num2,isCall::Bool=true) where {num1,num2<:Number}
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


function payoff(S::AbstractMatrix{num},euPayoff::EuropeanOption,spotData::equitySpotData,T1::num2=euPayoff.T) where{num,num2<:Number}
	r=spotData.r;
	T=euPayoff.T;
	iscall=euPayoff.isCall ? 1 : -1
	(Nsim,NStep)=size(S)
	NStep-=1;
	index1=round(Int,T/T1 * NStep)+1;
	ST=S[:,index1];
	K=euPayoff.K;
	f(ST::num__) where {num__<:Number}=max(iscall*(ST-K),0.0);
	payoff2=f.(ST);
	
	return payoff2*exp(-r*T);
end