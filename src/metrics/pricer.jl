	
function pricer(mcProcess::BaseProcess,spotData::equitySpotData,mcConfig::MonteCarloConfiguration,abstractPayoff::AbstractPayoff,mode1::MonteCarloMode=standard,parallelMode::BaseMode=SerialMode())
	Random.seed!(0)
	T=abstractPayoff.T;
	Nsim=mcConfig.Nsim;
	S=simulate(mcProcess,spotData,mcConfig,T,mode1,parallelMode)
	Payoff=payoff(S,abstractPayoff,spotData);
	Price=mean(Payoff);
	return Price;
end

function pricer(mcProcess::BaseProcess,spotData::equitySpotData,mcConfig::MonteCarloConfiguration,abstractPayoffs::Array{AbstractPayoff},mode1::MonteCarloMode=standard,parallelMode::BaseMode=SerialMode())
	Random.seed!(0)
	maxT=maximum([abstractPayoff.T for abstractPayoff in abstractPayoffs])
	Nsim=mcConfig.Nsim;
	S=simulate(mcProcess,spotData,mcConfig,maxT,mode1,parallelMode)
	Prices=[mean(payoff(S,abstractPayoff,spotData,maxT)) for abstractPayoff in abstractPayoffs  ]
	
	return Prices;
end
