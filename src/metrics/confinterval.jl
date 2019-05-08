using Distributions

"""
General Interface for Computation of confidence interval of price

		IC=confinter(mcProcess,spotData,mcBaseData,payoff_,,alpha,monteCarloMode=standard,parallelMode=SerialMode())
	
Where:\n
		mcProcess          = Process to be simulated.
		spotData  = Datas of the Spot.
		mcBaseData = Basic properties of MonteCarlo simulation
		payoff_ = Payoff(s) to be priced
		alpha [Optional, default to 99%] = confidence level
		monteCarloMode [Optional, default to standard]= standard or antitethic
		parallelMode  [Optional, default to SerialMode()] = SerialMode(), CudaMode(), AFMode()

		Price     = Price of the derivative

"""	
function confinter_macro(num1)
	@eval function confinter(mcProcess::$num1,spotData::equitySpotData,mcConfig::MonteCarloConfiguration,abstractPayoff::AbstractPayoff,alpha::Real=0.99,mode1::MonteCarloMode=standard,parallelMode::BaseMode=SerialMode())

		Random.seed!(0)
		T=abstractPayoff.T;
		Nsim=mcConfig.Nsim;
		S=simulate(mcProcess,spotData,mcConfig,T,mode1,parallelMode)
		Payoff=payoff(S,abstractPayoff,spotData);
		mean1=mean(Payoff);
		var1=var(Payoff);
		alpha_=1-alpha;
		dist1=Distributions.TDist(Nsim);
		tstar=quantile(dist1,1-alpha_/2.0)
		IC=(mean1-tstar*sqrt(var1)/Nsim,mean1+tstar*sqrt(var1)/Nsim)
		return IC;
		
	end
end

confinter_macro(BaseProcess)


function confinter_macro_array(num1)
	@eval function confinter(mcProcess::$num1,spotData::equitySpotData,mcConfig::MonteCarloConfiguration,abstractPayoffs::Array{AbstractPayoff},alpha::Real=0.99,mode1::MonteCarloMode=standard,parallelMode::BaseMode=SerialMode())

		Random.seed!(0)
		maxT=maximum([abstractPayoff.T for abstractPayoff in abstractPayoffs])
		Nsim=mcConfig.Nsim;
		S=simulate(mcProcess,spotData,mcConfig,maxT,mode1,parallelMode)
		Means=[mean(payoff(S,abstractPayoff,spotData,maxT)) for abstractPayoff in abstractPayoffs  ]
		Vars=[var(payoff(S,abstractPayoff,spotData,maxT)) for abstractPayoff in abstractPayoffs  ]
		alpha_=1-alpha;
		dist1=Distributions.TDist(Nsim);
		tstar=quantile(dist1,1-alpha_/2.0)
		IC=[(mean1-tstar*sqrt(var1)/Nsim,mean1+tstar*sqrt(var1)/Nsim) for (mean1,var1) in zip(Means,Vars)]
		return IC;
	end
end

confinter_macro_array(BaseProcess)