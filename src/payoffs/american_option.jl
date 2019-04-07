struct AmericanOption{num1,num2<:Number}<:AmericanPayoff
	T::num1
	K::num2
	isCall::Bool
	function AmericanOption(T::num1,K::num2,isCall::Bool=true) where {num1, num2 <: Number}
        if T <= 0.0
            error("Time to Maturity must be positive")
        elseif K <= 0.0
            error("Strike Price must be positive")
        else
            return new{num1,num2}(T,K,isCall)
        end
    end
end

export AmericanOption;

"""
Payoff computation from MonteCarlo paths

		Payoff=payoff(S,amPayoff,AmericanOption,)
	
Where:\n
		S           = Paths of the Underlying.
		amPayoff  = Datas of the Option.

		Payoff      = payoff of the option.
```
"""
function payoff(S::AbstractMatrix{numtype},amPayoff::AmericanOption,spotData::equitySpotData,T1::num2=amPayoff.T) where{numtype,num2<:Number}
	iscall=amPayoff.isCall ? 1 : -1
	K=amPayoff.K;
	T=amPayoff.T;
	isDualZero=typeof(zero(eltype(S))*K*0.0)
	phi_(Sti::numtype)::isDualZero=((Sti-K)*iscall>0.0) ? (Sti-K)*iscall : 0.0;
	(Nsim,NStep)=size(S)
	NStep-=1;
	index1=round(Int,T/T1 * NStep)+1;#round(Int,T/T1 * NStep)
	S1=view(S,:,1:index1)
	payoff1=payoff(collect(S1),spotData,phi_,T);
	
	return payoff1;
end
