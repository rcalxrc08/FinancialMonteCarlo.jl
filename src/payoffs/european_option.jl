
struct EuropeanOption<:AbstractEuropeanOptionData
	T::Float64
	K::Float64
	isCall::Bool
	function EuropeanOption(T::Float64,K::Float64,isCall::Bool=true)
        if T <= 0.0
            error("Time to Maturity must be positive")
        elseif K <= 0.0
            error("Strike Price must be positive")
        else
            return new(T,K,isCall)
        end
    end
end

export EuropeanOption;

"""
Payoff computation from MonteCarlo paths

		Payoff=payoff(S,euOptionData)
	
Where:\n
		S           = Paths of the Underlying.
		euOptionData  = Datas of the Option.

		Payoff      = payoff of the option.
```
"""
function payoff(S::Matrix{num},euOptionData::EuropeanOption,spotData::equitySpotData,T1::Float64=euOptionData.T) where{num<:Number}
	r=spotData.r;
	T=euOptionData.T;
	iscall=euOptionData.isCall?1:-1
	(Nsim,NStep)=size(S)
	NStep-=1;
	index1=round(Int,T/T1 * NStep)+1;
	S1=view(S,:,1:index1)
	ST=S1[:,end];
	K=euOptionData.K;
	f(ST::num)::num=(iscall*(ST-K)>0.0)?iscall*(ST-K):0.0;
	payoff2=f.(ST);
	
	return payoff2*exp(-r*T);
end