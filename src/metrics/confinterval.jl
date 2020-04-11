using Distributions

"""
General Interface for Computation of confidence interval of price

		IC=confinter(mcProcess,rfCurve,mcBaseData,payoff_,alpha)
	
Where:\n
		mcProcess          = Process to be simulated.
		rfCurve  = Zero Rate Data.
		mcBaseData = Basic properties of MonteCarlo simulation
		payoff_ = Payoff(s) to be priced
		alpha [Optional, default to 99%] = confidence level
		

		Price     = Price of the derivative

"""	
function confinter(mcProcess::BaseProcess,rfCurve::AbstractZeroRateCurve,mcConfig::MonteCarloConfiguration,abstractPayoff::AbstractPayoff,alpha::Real=0.99)

	set_seed(mcConfig)
	T=maturity(abstractPayoff);
	Nsim=mcConfig.Nsim;
	S=simulate(mcProcess,rfCurve,mcConfig,T)
	Payoff=payoff(S,abstractPayoff,rfCurve);
	mean1=mean(Payoff);
	var1=var(Payoff);
	alpha_=1-alpha;
	dist1=Distributions.TDist(Nsim);
	tstar=quantile(dist1,1-alpha_/2.0)
	IC=(mean1-tstar*sqrt(var1)/Nsim,mean1+tstar*sqrt(var1)/Nsim)
	return IC;
	
end


function confinter(mcProcess::BaseProcess,rfCurve::AbstractZeroRateCurve,mcConfig::MonteCarloConfiguration,abstractPayoffs::Array{abstractPayoff_},alpha::Real=0.99) where {abstractPayoff_ <: AbstractPayoff}

	set_seed(mcConfig)
	maxT=maximum([maturity(abstractPayoff) for abstractPayoff in abstractPayoffs])
	Nsim=mcConfig.Nsim;
	S=simulate(mcProcess,rfCurve,mcConfig,maxT)
	Means=[mean(payoff(S,abstractPayoff,rfCurve,maxT)) for abstractPayoff in abstractPayoffs  ]
	Vars=[var(payoff(S,abstractPayoff,rfCurve,maxT)) for abstractPayoff in abstractPayoffs  ]
	alpha_=1-alpha;
	dist1=Distributions.TDist(Nsim);
	tstar=quantile(dist1,1-alpha_/2.0)
	IC=[(mean1-tstar*sqrt(var1)/Nsim,mean1+tstar*sqrt(var1)/Nsim) for (mean1,var1) in zip(Means,Vars)]
	return IC;
end
