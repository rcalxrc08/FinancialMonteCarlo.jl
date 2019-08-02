
function pricer_macro(model_type)
	@eval begin
		"""
		General Interface for Pricing

				Price=pricer(mcProcess,spotData,mcBaseData,payoff_,parallelMode=SerialMode())
			
		Where:\n
				mcProcess          = Process to be simulated.
				spotData  = Datas of the Spot.
				mcBaseData = Basic properties of MonteCarlo simulation
				payoff_ = Payoff(s) to be priced
				parallelMode  [Optional, default to SerialMode()] = SerialMode(), CudaMode(), AFMode()

				Price     = Price of the derivative

		"""	
		function pricer(mcProcess::$model_type,spotData::equitySpotData,mcConfig::MonteCarloConfiguration,abstractPayoff::AbstractPayoff,parallelMode::BaseMode=SerialMode())

			Random.seed!(0)
			T=abstractPayoff.T;
			S=simulate(mcProcess,spotData,mcConfig,T,parallelMode)
			Payoff=payoff(S,abstractPayoff,spotData);
			Price=mean(Payoff);
			return Price;
		end
	end
end

pricer_macro(BaseProcess)

function pricer_macro_array(model_type)
	@eval function pricer(mcProcess::$model_type,spotData::equitySpotData,mcConfig::MonteCarloConfiguration,abstractPayoffs::Array{AbstractPayoff},parallelMode::BaseMode=SerialMode())
		Random.seed!(0)
		maxT=maximum([abstractPayoff.T for abstractPayoff in abstractPayoffs])
		S=simulate(mcProcess,spotData,mcConfig,maxT,parallelMode)
		Prices=[mean(payoff(S,abstractPayoff,spotData,maxT)) for abstractPayoff in abstractPayoffs  ]
		
		return Prices;
	end
end

pricer_macro_array(BaseProcess)