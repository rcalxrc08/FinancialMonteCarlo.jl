type Forward<:EuropeanPayoff end

struct ForwardData<:AbstractEuropeanOptionData
	T::Float64
end

export Forward,ForwardData;

function payoff(S::Matrix{num},optionData::ForwardData,Payoff::Forward,isCall::Bool=true) where{num<:Number}
	
	ST=S[1:end,end];
	
	return ST;
end