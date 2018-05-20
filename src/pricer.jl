include("utils.jl")
include("payoff.jl")
include("models.jl")
	
function pricer(mcProcess::AbstractMonteCarloProcess,spotData::equitySpotData,mcBaseData::MonteCarloBaseData,optionData::OptionData,payoff1::PayoffMC,isCall::Bool=true,mode1::MonteCarloMode=standard)
	srand(0)
	T=optionData.T;
	Nsim=mcBaseData.Nsim;
	S=simulate(mcProcess,spotData,mcBaseData,T,mode1)
	Payoff=payoff(S,optionData,spotData,payoff1,isCall);
	Price=mean(Payoff);
	return Price;
end

function pricer(mcProcess::MonteCarloProblem,spotData::equitySpotData,mcBaseData::MonteCarloBaseData,optionData::OptionData,payoff1::PayoffMC,isCall::Bool=true,mode1::MonteCarloMode=standard)
	srand(0)
	T=optionData.T;
	Nsim=mcBaseData.Nsim;
	S=simulate(mcProcess,spotData,mcBaseData,T,mode1)
	Payoff=payoff(S,optionData,spotData,payoff1,isCall);
	Price=mean(Payoff);
	return Price;
end


function pricer(mcProcess::AbstractMonteCarloProcess,spotData::equitySpotData,mcBaseData::MonteCarloBaseData,optionDatas::Array{OptionData},payoffs1::Array{PayoffMC},isCalls::BitArray{1}=trues(length(optionDatas)),mode1::MonteCarloMode=standard)
	srand(0)
	maxT=maximum([optionData.T for optionData in optionDatas])
	r=spotData.r;
	Nsim=mcBaseData.Nsim;
	S=simulate(mcProcess,spotData,mcBaseData,maxT,mode1)
	Prices=[mean(payoff(S,optionData,spotData,payoff1,isCall,maxT)) for (payoff1,optionData,isCall) in zip(payoffs1,optionDatas,isCalls)  ]
	
	return Prices;
end

function pricer(mcProcess::MonteCarloProblem,spotData::equitySpotData,mcBaseData::MonteCarloBaseData,optionDatas::Array{OptionData},payoffs1::Array{PayoffMC},isCalls::BitArray{1}=trues(length(optionDatas)),mode1::MonteCarloMode=standard)
	srand(0)
	maxT=maximum([optionData.T for optionData in optionDatas])
	Nsim=mcBaseData.Nsim;
	S=simulate(mcProcess,spotData,mcBaseData,maxT,mode1)
	Prices=[mean(payoff(S,optionData,spotData,payoff1,isCall,maxT)) for (payoff1,optionData,isCall) in zip(payoffs1,optionDatas,isCalls)  ]
	
	return Prices;
end
