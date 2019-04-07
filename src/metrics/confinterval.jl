using Distributions
	
function confinter(mcProcess::BaseProcess,spotData::equitySpotData,mcConfig::MonteCarloConfiguration,abstractPayoff::AbstractPayoff,mode1::MonteCarloMode=standard,alpha::Real=0.99,parallelMode::BaseMode=SerialMode())
	Random.seed!(0)
	T=abstractPayoff.T;
	Nsim=mcConfig.Nsim;
	S=simulate(mcProcess,spotData,mcConfig,T,mode1)
	Payoff=payoff(S,abstractPayoff,spotData);
	mean1=mean(Payoff);
	var1=var(Payoff);
	dist1=Distributions.TDist(Nsim);
	tstar=quantile(dist1,1-alpha/2.0)
	IC=(mean1-tstar*sqrt(var1)/Nsim,mean1+tstar*sqrt(var1)/Nsim)
	return IC;
end

function confinter(mcProcess::BaseProcess,spotData::equitySpotData,mcConfig::MonteCarloConfiguration,abstractPayoffs::Array{AbstractPayoff},mode1::MonteCarloMode=standard,alpha::Real=0.99,parallelMode::BaseMode=SerialMode())
	Random.seed!(0)
	maxT=maximum([abstractPayoff.T for abstractPayoff in abstractPayoffs])
	Nsim=mcConfig.Nsim;
	S=simulate(mcProcess,spotData,mcConfig,maxT,mode1)
	Means=[mean(payoff(S,abstractPayoff,spotData,maxT)) for abstractPayoff in abstractPayoffs  ]
	Vars=[var(payoff(S,abstractPayoff,spotData,maxT)) for abstractPayoff in abstractPayoffs  ]
	dist1=Distributions.TDist(Nsim);
	tstar=quantile(dist1,1-alpha/2.0)
	IC=[(mean1-tstar*sqrt(var1)/Nsim,mean1+tstar*sqrt(var1)/Nsim) for (mean1,var1) in zip(Means,Vars)]
	return IC;
end
