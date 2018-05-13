include("utils.jl")
include("payoff.jl")
include("models.jl")
	
function pricer(mcProcess::AbstractMonteCarloProcess,spotData::equitySpotData,mcBaseData::MonteCarloBaseData,optionData::OptionData,payoff1::PayoffMC,isCall::Bool=true,mode1::MonteCarloMode=standard)
	srand(0)
	T=optionData.T;
	r=spotData.r;
	Nsim=mcBaseData.Nsim;
	S=simulate(mcProcess,spotData,mcBaseData,T)
	Payoff=payoff(S,optionData,spotData,payoff1,isCall);
	Price=mean(Payoff);
	return Price;
end

