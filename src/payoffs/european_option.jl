
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

"""
Payoff computation from MonteCarlo paths

		Payoff=payoff(S,euPayoff)
	
Where:\n
		S           = Paths of the Underlying.
		euPayoff  = Datas of the Option.

		Payoff      = payoff of the option.
```
"""
function payoff(S::AbstractMatrix{num},euPayoff::EuropeanOption,spotData::equitySpotData,T1::num2=euPayoff.T) where{num,num2<:Number}
	r=spotData.r;
	T=euPayoff.T;
	iscall=euPayoff.isCall ? 1 : -1
	(Nsim,NStep)=size(S)
	NStep-=1;
	index1=round(Int,T/T1 * NStep)+1;
	S1=view(S,:,1:index1)
	ST=S1[:,end];
	K=euPayoff.K;
	f(ST::num)::num=(iscall*(ST-K)>0.0)?iscall*(ST-K) : 0.0;
	payoff2=f.(ST);
	
	return payoff2*exp(-r*T);
end