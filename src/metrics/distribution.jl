
function distribution_macro(model_type)
	@eval begin
		"""
		General Interface for Pricing

				Price=pricer(mcProcess,spotData,mcBaseData,payoff_)
			
		Where:\n
				mcProcess          = Process to be simulated.
				spotData  = Datas of the Spot.
				mcBaseData = Basic properties of MonteCarlo simulation
				payoff_ = Payoff(s) to be priced
				

				Price     = Price of the derivative

		"""	
		function distribution(mcProcess::$model_type,spotData::ZeroRateCurve,mcConfig::MonteCarloConfiguration)

			set_seed(mcConfig)
			T=maturity(abstractPayoff);
			S=simulate(mcProcess,spotData,mcConfig,T)
	
			return S[:,end];
		end
	end
end

distribution_macro(BaseProcess)