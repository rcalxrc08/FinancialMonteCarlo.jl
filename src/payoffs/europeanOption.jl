type EuropeanOption<:PayoffMC end

struct EUOptionData<:OptionData
	T::Float64
	K::Float64
end

function payoff(S::Matrix{num},optionData::EUOptionData,Payoff::EuropeanOption,isCall::Bool=true) where{num<:Number}
	iscall=isCall?1:-1
	ST=S[1:end,end];
	K=optionData.K;
	f(ST::num)::num=(iscall*(ST-K)>0.0)?iscall*(ST-K):0.0;
	payoff2=f.(ST);
	
	return payoff2;
end