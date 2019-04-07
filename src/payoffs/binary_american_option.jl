struct BinaryAmericanOption{num1,num2<:Number}<:AmericanPayoff
	T::num1
	K::num2
	isCall::Bool
	function BinaryAmericanOption(T::num1,K::num2,isCall::Bool=true) where {num1,num2<:Number}
        if T <= 0.0
            error("Time to Maturity must be positive")
        elseif K <= 0.0
            error("Strike Price must be positive")
        else
            return new{num1,num2}(T,K,isCall)
        end
    end
end

export BinaryAmericanOption;

"""
Payoff computation from MonteCarlo paths

		Payoff=payoff(S,amPayoff,BinaryAmericanOption,)
	
Where:\n
		S           = Paths of the Underlying.
		amPayoff  = Datas of the Option.
		BinaryAmericanOption = Type of the Option

		Payoff      = payoff of the option.
```
"""
function payoff(S::AbstractMatrix{num},amPayoff::BinaryAmericanOption,spotData::equitySpotData,T1::num2=amPayoff.T) where{num,num2<:Number}
	iscall=amPayoff.isCall ? 1 : -1
	K=amPayoff.K;
	T=amPayoff.T;
	(Nsim,NStep)=size(S)
	NStep-=1;
	index1=round(Int,T/T1 * NStep)+1;
	S1=view(S,:,1:index1)
	isDualZero=typeof(S[1,1]*K*0.0)
	phi(Sti::num)::isDualZero where {num<:Number}=((Sti-K)*iscall>0.0) ? 1.0 : 0.0;
	
	payoff1=payoff(collect(S1),spotData,phi,T);
	
	return payoff1;
end
