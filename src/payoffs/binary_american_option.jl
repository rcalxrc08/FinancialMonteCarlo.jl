"""
Class for Dispatching Forward Payoff

		forward=Forward(T::num) where {num<:Number}
	
Where:\n
		T	=	Time to maturity of the Forward.
"""
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


function payoff(S::AbstractMatrix{num},amPayoff::BinaryAmericanOption,spotData::equitySpotData,T1::num2=amPayoff.T) where{num,num2<:Number}
	iscall=amPayoff.isCall ? 1 : -1
	K=amPayoff.K;
	T=amPayoff.T;
	(Nsim,NStep)=size(S)
	NStep-=1;
	index1=round(Int,T/T1 * NStep)+1;
	S1=view(S,:,1:index1)
	isDualZero=typeof(S[1,1]*K*0.0)
	phi_(Sti::num)::isDualZero=((Sti-K)*iscall>0.0) ? 1.0 : 0.0;
	
	payoff1=payoff(collect(S1),spotData,phi_,T);
	
	return payoff1;
end
