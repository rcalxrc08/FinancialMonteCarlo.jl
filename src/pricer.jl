include("payoff.jl")
include("struct_utils.jl")
include("models.jl")
	
function pricer(mcProcess::AbstractMonteCarloProcess,spotData::equitySpotData,mcBaseData::MonteCarloBaseData,optionData::OptionData,payoff1::PayoffMC,isCall::Bool=true)
	srand(0)
	T=optionData.T;
	r=spotData.r;
	Nsim=mcBaseData.Nsim;
	S=simulate(mcProcess,spotData,mcBaseData,T)
	Payoff=payoff(S,optionData,payoff1,isCall);
	Price=mean(Payoff)*exp(-r*T);
	return Price;
end

