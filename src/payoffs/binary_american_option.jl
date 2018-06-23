struct BinaryAmericanOption<:AbstractEuropeanOptionData
	T::Float64
	K::Float64
	function BinaryAmericanOption(T::Float64,K::Float64)
        if T <= 0.0
            error("Time to Maturity must be positive")
        elseif K <= 0.0
            error("Strike Price must be positive")
        else
            return new(T,K)
        end
    end
end

"""
Payoff computation from MonteCarlo paths

		Payoff=payoff(S,amOptionData,BinaryAmericanOption,isCall=true)
	
Where:\n
		S           = Paths of the Underlying.
		amOptionData  = Datas of the Option.
		BinaryAmericanOption = Type of the Option
		isCall = true for Call Options, false for price Options.

		Payoff      = payoff of the option.
```
"""
function payoff(S::Matrix{num},amOptionData::BinaryAmericanOption,spotData::equitySpotData,isCall::Bool=true,T1::Float64=amOptionData.T) where{num<:Number}
	iscall=isCall?1:-1
	K=amOptionData.K;
	T=amOptionData.T;
	(Nsim,NStep)=size(S)
	NStep-=1;
	index1=round(Int,T/T1 * NStep)+1;
	S1=view(S,:,1:index1)
	phi(Sti::Number)::Number=((Sti-K)*iscall>0.0)?1.0:0.0;
	
	payoff1=payoff(collect(S1),spotData,phi,T);
	
	return payoff1;
end
