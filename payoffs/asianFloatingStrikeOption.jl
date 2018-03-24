type AsianFloatingStrikeOption<:PayoffMC end

struct AsianFloatingStrikeOptionData<:OptionData
	T::Float64
end

function payoff(S::Matrix{num},optionData::AsianFloatingStrikeOptionData,Payoff::AsianFloatingStrikeOption,isCall::Bool=true) where{num<:Number}
	iscall=isCall?1:-1
	NsimTmp=length(S[1:end,end]);
	@inbounds f(S::Array{num})::num=(iscall*(S[end]-mean(S))>0.0)?(iscall*(S[end]-mean(S))):0.0;		
	@inbounds payoff2=[f(S[i,1:end]) for i in 1:NsimTmp];
	
	return payoff2;
end