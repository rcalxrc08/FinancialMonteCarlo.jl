struct AmericanOption<:AbstractEuropeanOptionData
	T::Float64
	K::Float64
	isCall::Bool
	function AmericanOption(T::Float64,K::Float64,isCall::Bool=true)
        if T <= 0.0
            error("Time to Maturity must be positive")
        elseif K <= 0.0
            error("Strike Price must be positive")
        else
            return new(T,K,isCall)
        end
    end
end

export AmericanOption;

"""
Payoff computation from MonteCarlo paths

		Payoff=payoff(S,amOptionData,AmericanOption,)
	
Where:\n
		S           = Paths of the Underlying.
		amOptionData  = Datas of the Option.

		Payoff      = payoff of the option.
```
"""
function payoff(S::Matrix{num},amOptionData::AmericanOption,spotData::equitySpotData,T1::Float64=amOptionData.T) where{num<:Number}
	iscall=amOptionData.isCall?1:-1
	K=amOptionData.K;
	T=amOptionData.T;
	phi(Sti::Number)::Number=((Sti-K)*iscall>0.0)?(Sti-K)*iscall:0.0;
	(Nsim,NStep)=size(S)
	NStep-=1;
	index1=round(Int,T/T1 * NStep)+1;#round(Int,T/T1 * NStep)
	S1=view(S,:,1:index1)
	payoff1=payoff(collect(S1),spotData,phi,T);
	
	return payoff1;
end
