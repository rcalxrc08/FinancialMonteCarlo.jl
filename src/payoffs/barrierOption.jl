type BarrierOption<:PayoffMC end

struct BarrierOptionData<:OptionData
	T::Float64
	K::Float64
	D::Float64
end

export BarrierOption,BarrierOptionData;

function payoff(S::Matrix{num},optionData::BarrierOptionData,Payoff::BarrierOption,isCall::Bool=true) where{num<:Number}
	iscall=isCall?1:-1
	NsimTmp=length(S[1:end,end]);
	K=optionData.K;
	D=optionData.D;
	@inbounds f(S::Array{num})::num=(iscall*(S[end]-K)>0.0)&&(minimum(S)>D)?iscall*(S[end]-K):0.0;		
	@inbounds payoff2=[f(S[i,1:end]) for i in 1:NsimTmp];
	
	return payoff2;
end