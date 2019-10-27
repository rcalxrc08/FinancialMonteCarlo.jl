
function pricer_macro(model_type)
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
		function pricer(mcProcess::$model_type,spotData::equitySpotData,mcConfig::MonteCarloConfiguration,abstractPayoff::AbstractPayoff)

			set_seed(mcConfig)
			T=abstractPayoff.T;
			S=simulate(mcProcess,spotData,mcConfig,T)
			Payoff=payoff(S,abstractPayoff,spotData);
			Price=mean(Payoff);
			return Price;
		end
	end
end

pricer_macro(BaseProcess)

function pricer_macro_array(model_type)
	@eval function pricer(mcProcess::$model_type,spotData::equitySpotData,mcConfig::MonteCarloConfiguration,abstractPayoffs::Array{abstractPayoff}) where {abstractPayoff <: AbstractPayoff}
		set_seed(mcConfig)
		maxT=maximum([abstractPayoff.T for abstractPayoff in abstractPayoffs])
		S=simulate(mcProcess,spotData,mcConfig,maxT)
		Prices=[mean(payoff(S,abstractPayoff,spotData,maxT)) for abstractPayoff in abstractPayoffs  ]
		
		return Prices;
	end
end

function pricer_macro_dict(model_type)
	@eval function pricer(mcProcess::$model_type,spotData::equitySpotData,mcConfig::MonteCarloConfiguration,dict_::Dict{FinancialMonteCarlo.AbstractPayoff,Number})
		set_seed(mcConfig)
		abstractPayoffs=keys(dict_);
		weights_=values(dict_);
		maxT=maximum([abstractPayoff.T for abstractPayoff in abstractPayoffs])
		S=simulate(mcProcess,spotData,mcConfig,maxT)
		price=sum(weight_*mean(payoff(S,abstractPayoff,spotData,maxT)) for (abstractPayoff,weight_) in dict_);
		
		return price;
	end
end

pricer_macro_array(BaseProcess)
pricer_macro_dict(BaseProcess)