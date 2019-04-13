"""
Struct for Standard American Option

		AmericanOption{num1,num2<:Number}<:AmericanPayoff
	
Fields:\n
		T           = Time to Maturity of the Options.
		K  = Strike of the Option.
		isCall  = true for call, false for put.
```
"""
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


function payoff(S::AbstractMatrix{numtype},amPayoff::AmericanOption,spotData::equitySpotData,T1::num2=amPayoff.T) where{numtype,num2<:Number}
	iscall=amPayoff.isCall ? 1 : -1
	K=amPayoff.K;
	T=amPayoff.T;
	isDualZero=typeof(zero(eltype(S))*K)
	phi_(Sti::numtype_) where {numtype_<:Number}=((Sti-K)*iscall>0.0) ? (Sti-K)*iscall : zero(isDualZero);
	(Nsim,NStep)=size(S)
	NStep-=1;
	index1=round(Int,T/T1 * NStep)+1;#round(Int,T/T1 * NStep)
	S1=view(S,:,1:index1)
	payoff1=payoff(collect(S1),spotData,phi_,T);
	
	return payoff1;
end
