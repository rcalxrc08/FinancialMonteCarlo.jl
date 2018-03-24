type Forward<:PayoffMC end

struct ForwardData<:OptionData
	T::Float64
end

export Forward,ForwardData;

function payoff(S::Matrix{num},optionData::ForwardData,Payoff::Forward,isCall::Bool=true) where{num<:Number}
	ST=S[1:end,end];
	f(ST::num)::num=ST;
	payoff2=f.(ST);
	
	return payoff2;
end