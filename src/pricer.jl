include("utils.jl")
include("payoff.jl")
include("models.jl")
	
function pricer(mcProcess::AbstractMonteCarloProcess,spotData::equitySpotData,mcBaseData::MonteCarloBaseData,optionData::OptionData,payoff1::PayoffMC,isCall::Bool=true,mode1::MonteCarloMode=standard)
	srand(0)
	T=optionData.T;
	r=spotData.r;
	Nsim=mcBaseData.Nsim;
	S=simulate(mcProcess,spotData,mcBaseData,T,mode1)
	Payoff=payoff(S,optionData,spotData,payoff1,isCall);
	Price=mean(Payoff);
	return Price;
end

function pricer(mcProcess::MonteCarloProblem,spotData::equitySpotData,mcBaseData::MonteCarloBaseData,optionData::OptionData,payoff1::PayoffMC,isCall::Bool=true,mode1::MonteCarloMode=standard)
	srand(0)
	T=optionData.T;
	r=spotData.r;
	Nsim=mcBaseData.Nsim;
	S=simulate(mcProcess,spotData,mcBaseData,T,mode1)
	Payoff=payoff(S,optionData,spotData,payoff1,isCall);
	Price=mean(Payoff);
	return Price;
end


function pricer(mcProcess::MonteCarloProblem,spotData::equitySpotData,mcBaseData::MonteCarloBaseData,optionDatas::Array{OptionData},payoffs1::Array{PayoffMC},isCalls::Array{Bool},mode1::MonteCarloMode=standard)
	srand(0)
	maxT=-1.0;
	for optionData in optionDatas
		maxT=maxT<T?T:maxT;
	end
	T=optionDatas.T;
	r=spotData.r;
	Nsim=mcBaseData.Nsim;
	S=simulate(mcProcess,spotData,mcBaseData,maxT,mode1)
	Prices=[mean(payoff(S,optionData,spotData,payoff1,isCall,maxT)) for (payoff1,optionData,isCall) in zip(payoffs1,optionDatas,isCalls)  ]
	
	return Prices;
end

