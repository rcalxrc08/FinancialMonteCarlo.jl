type Forward<:EuropeanPayoff end

struct ForwardData<:AbstractEuropeanOptionData
	T::Float64
end

export Forward,ForwardData;

function payoff(S::Matrix{num},optionData::ForwardData,spotData::equitySpotData,Payoff::Forward,isCall::Bool=true,T1::Float64=T) where{num<:Number}
	r=spotData.r;
	T=optionData.T;
	NStep=length(S[1,1:end])-1
	dt=T1/NStep;
	index1=round(Int,T/dt);
	index1=index1>NStep?Nstep+1:index+1;
	ST=S[1:end,index1+1];
	
	return ST*exp(-r*T);
end