	
function pricer(mcProcess::BaseProcess,spotData::equitySpotData,mcConfig::MonteCarloConfiguration,abstractPayoff::AbstractPayoff,mode1::MonteCarloMode=standard)
	srand(0)
	T=abstractPayoff.T;
	Nsim=mcConfig.Nsim;
	S=simulate(mcProcess,spotData,mcConfig,T,mode1)
	Payoff=payoff(S,abstractPayoff,spotData);
	Price=mean(Payoff);
	return Price;
end

function pricer(mcProcess::BaseProcess,spotData::equitySpotData,mcConfig::MonteCarloConfiguration,abstractPayoffs::Array{AbstractPayoff},mode1::MonteCarloMode=standard)
	srand(0)
	maxT=maximum([abstractPayoff.T for abstractPayoff in abstractPayoffs])
	Nsim=mcConfig.Nsim;
	S=simulate(mcProcess,spotData,mcConfig,maxT,mode1)
	Prices=[mean(payoff(S,abstractPayoff,spotData,maxT)) for abstractPayoff in abstractPayoffs  ]
	
	return Prices;
end
