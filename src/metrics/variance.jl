
function variance_macro(model_type)
	@eval begin
		"""
		General Interface for Computation of variance interval of price

				variance_=variance(mcProcess,spotData,mcBaseData,payoff_,)
			
		Where:\n
				mcProcess          = Process to be simulated.
				spotData  = Datas of the Spot.
				mcBaseData = Basic properties of MonteCarlo simulation
				payoff_ = Payoff(s) to be priced
				

				variance_     = variance of the payoff of the derivative

		"""
		function variance(mcProcess::$model_type,spotData::equitySpotData,mcConfig::MonteCarloConfiguration,abstractPayoff::AbstractPayoff)
			set_seed(mcConfig)
			T=maturity(abstractPayoff);
			S=simulate(mcProcess,spotData,mcConfig,T)
			Payoff=payoff(S,abstractPayoff,spotData);
			variance_=var(Payoff);
			return variance_;
		end
	end
end

variance_macro(BaseProcess)


function variance_macro_array(model_type)
	@eval function variance(mcProcess::$model_type,spotData::equitySpotData,mcConfig::MonteCarloConfiguration,abstractPayoffs::Array{AbstractPayoff})
		set_seed(mcConfig)
		maxT=maximum([maturity(abstractPayoff) for abstractPayoff in abstractPayoffs])
		S=simulate(mcProcess,spotData,mcConfig,maxT)
		variance_=[var(payoff(S,abstractPayoff,spotData,maxT)) for abstractPayoff in abstractPayoffs  ]
		
		return variance_;
	end
end

variance_macro_array(BaseProcess)