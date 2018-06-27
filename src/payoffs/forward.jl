
struct Forward<:EuropeanPayoff
	T::Float64
	function Forward(T::Float64)
        if T <= 0.0
            error("Time to Maturity must be positive")
        else
            return new(T)
        end
    end
end

export Forward;

function payoff(S::AbstractMatrix{num},optionData::Forward,spotData::equitySpotData,T1::Float64=optionData.T) where{num<:Number}
	r=spotData.r;
	T=optionData.T;
	(Nsim,NStep)=size(S)
	NStep-=1;
	index1=round(Int,T/T1 * NStep)+1;
	S1=view(S,:,1:index1)
	ST=S1[:,end];
	
	return ST*exp(-r*T);
end