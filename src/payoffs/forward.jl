type Forward<:EuropeanPayoff end

struct ForwardData<:AbstractEuropeanOptionData
	T::Float64
end

export Forward,ForwardData;

function payoff(S::Matrix{num},optionData::ForwardData,spotData::equitySpotData,Payoff::Forward,isCall::Bool=true) where{num<:Number}
	r=spotData.r;
	T=optionData.T;
	ST=S[1:end,end];
	
	return ST*exp(-r*T);
end