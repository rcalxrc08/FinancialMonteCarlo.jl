type AsianFixedStrikeOption<:PayoffMC end

struct AsianFixedStrikeOptionData<:OptionData
	T::Float64
	K::Float64
end

export AsianFixedStrikeOption,AsianFixedStrikeOptionData;

function payoff(S::Matrix{num},optionData::AsianFixedStrikeOptionData,Payoff::AsianFixedStrikeOption,isCall::Bool=true) where{num<:Number}
	iscall=isCall?1:-1
	K=optionData.K;
	NsimTmp=length(S[1:end,end]);
	@inbounds f(S::Array{num})::num=(iscall*(mean(S)-K)>0.0)?(iscall*(mean(S)-K)):0.0;		
	@inbounds payoff2=[f(S[i,1:end]) for i in 1:NsimTmp];
	
	return payoff2;
end