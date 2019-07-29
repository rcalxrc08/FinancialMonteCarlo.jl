
function variance_macro(num1)
	@eval begin
		"""
		General Interface for Computation of variance interval of price

				variance_=variance(mcProcess,spotData,mcBaseData,payoff_,,monteCarloMode=standard,parallelMode=SerialMode())
			
		Where:\n
				mcProcess          = Process to be simulated.
				spotData  = Datas of the Spot.
				mcBaseData = Basic properties of MonteCarlo simulation
				payoff_ = Payoff(s) to be priced
				monteCarloMode [Optional, default to standard]= standard or antitethic
				parallelMode  [Optional, default to SerialMode()] = SerialMode(), CudaMode(), AFMode()

				variance_     = variance of the payoff of the derivative

		"""
		function variance(mcProcess::$num1,spotData::equitySpotData,mcConfig::MonteCarloConfiguration,abstractPayoff::AbstractPayoff,parallelMode::BaseMode=SerialMode())
			Random.seed!(0)
			T=abstractPayoff.T;
			S=simulate(mcProcess,spotData,mcConfig,T,parallelMode)
			Payoff=payoff(S,abstractPayoff,spotData);
			variance_=var(Payoff);
			return variance_;
		end
	end
end

variance_macro(BaseProcess)


function variance_macro_array(num1)
	@eval function variance(mcProcess::$num1,spotData::equitySpotData,mcConfig::MonteCarloConfiguration,abstractPayoffs::Array{AbstractPayoff},parallelMode::BaseMode=SerialMode())
		Random.seed!(0)
		maxT=maximum([abstractPayoff.T for abstractPayoff in abstractPayoffs])
		S=simulate(mcProcess,spotData,mcConfig,maxT,parallelMode)
		variance_=[var(payoff(S,abstractPayoff,spotData,maxT)) for abstractPayoff in abstractPayoffs  ]
		
		return variance_;
	end
end

variance_macro_array(BaseProcess)