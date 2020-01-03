
function pricer_macro(model_type)
	@eval begin
		"""
		General Interface for Pricing

				Price=pricer(mcProcess,rfCurve,mcBaseData,payoff_)
			
		Where:\n
				mcProcess          = Process to be simulated.
				rfCurve  = Datas of the Spot.
				mcBaseData = Basic properties of MonteCarlo simulation
				payoff_ = Payoff(s) to be priced
				

				Price     = Price of the derivative

		"""	
		function pricer(mcProcess::$model_type,rfCurve::AbstractZeroRateCurve,mcConfig::MonteCarloConfiguration,abstractPayoff::AbstractPayoff)

			set_seed(mcConfig)
			T=maturity(abstractPayoff);
			S=simulate(mcProcess,rfCurve,mcConfig,T)
			Payoff=payoff(S,abstractPayoff,rfCurve);
			Price=mean(Payoff);
			return Price;
		end
	end
end

pricer_macro(BaseProcess)

function pricer_macro_array(model_type)
	@eval function pricer(mcProcess::$model_type,rfCurve::AbstractZeroRateCurve,mcConfig::MonteCarloConfiguration,abstractPayoffs::Array{abstractPayoff_}) where {abstractPayoff_ <: AbstractPayoff}
		set_seed(mcConfig)
		maxT=maximum([maturity(abstractPayoff) for abstractPayoff in abstractPayoffs])
		S=simulate(mcProcess,rfCurve,mcConfig,maxT)
		Prices=[mean(payoff(S,abstractPayoff,rfCurve,maxT)) for abstractPayoff in abstractPayoffs  ]
		
		return Prices;
	end
end

function pricer_macro_dict(model_type)
	@eval function pricer(mcProcess::$model_type,rfCurve::AbstractZeroRateCurve,mcConfig::MonteCarloConfiguration,dict_::Dict{FinancialMonteCarlo.AbstractPayoff,Number})
		set_seed(mcConfig)
		abstractPayoffs=keys(dict_);
		maxT=maximum([maturity(abstractPayoff) for abstractPayoff in abstractPayoffs])
		S=simulate(mcProcess,rfCurve,mcConfig,maxT)
		price=sum(weight_*mean(payoff(S,abstractPayoff,rfCurve,maxT)) for (abstractPayoff,weight_) in dict_);
		
		return price;
	end
end

function pricer(mcProcess::Dict{String,FinancialMonteCarlo.AbstractMonteCarloProcess},rfCurve::AbstractZeroRateCurve,mcConfig::MonteCarloConfiguration,dict_::Dict{String,Dict{FinancialMonteCarlo.AbstractPayoff,Number}})
		set_seed(mcConfig)
		underlyings_payoff=keys(dict_)
		price=0.0;
		for under_ in underlyings_payoff
			options=dict_[under_]
			model=mcProcess[under_]
			price=price+pricer(model,rfCurve,mcConfig,options);
		end
		
		return price;
	end

pricer_macro_array(BaseProcess)
pricer_macro_dict(BaseProcess)