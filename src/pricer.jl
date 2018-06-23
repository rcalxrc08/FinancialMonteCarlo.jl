include("utils.jl")
include("payoff.jl")
include("models.jl")
	
function pricer(mcProcess::BaseProcess,spotData::equitySpotData,mcBaseData::MonteCarloBaseData,optionData::OptionData,mode1::MonteCarloMode=standard)
	srand(0)
	T=optionData.T;
	Nsim=mcBaseData.Nsim;
	S=simulate(mcProcess,spotData,mcBaseData,T,mode1)
	Payoff=payoff(S,optionData,spotData);
	Price=mean(Payoff);
	return Price;
end

function pricer(mcProcess::BaseProcess,spotData::equitySpotData,mcBaseData::MonteCarloBaseData,optionDatas::Array{OptionData},mode1::MonteCarloMode=standard)
	srand(0)
	maxT=maximum([optionData.T for optionData in optionDatas])
	Nsim=mcBaseData.Nsim;
	S=simulate(mcProcess,spotData,mcBaseData,maxT,mode1)
	Prices=[mean(payoff(S,optionData,spotData,maxT)) for optionData in optionDatas  ]
	
	return Prices;
end
