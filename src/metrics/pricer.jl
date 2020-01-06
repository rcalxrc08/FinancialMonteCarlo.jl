include("utils.jl")
function pricer_macro(model_type)
	@eval begin
		"""
		Low Level Interface for Pricing

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

#####Pricer for multivariate
function pricer(mcProcess::VectorialMonteCarloProcess,rfCurve::AbstractZeroRateCurve,mcConfig::MonteCarloConfiguration,dict_::Array{Dict{FinancialMonteCarlo.AbstractPayoff,Number}})
	set_seed(mcConfig)
	N_=length(mcProcess.models);
	idx_=compute_indices(N_);
	## Filter dict_ for undef
	IND_=collect(1:length(dict_));
	filter!(i->isassigned(dict_,i),IND_)
	dict_cl=dict_[IND_];
	
	maxT=maximum([ maximum(FinancialMonteCarlo.maturity.(collect(keys(ar_el)))) for ar_el in dict_cl])
	S=simulate(mcProcess,rfCurve,mcConfig,maxT)
	#price=0.0;
	#for i in IND_
	#	idx_loc=idx_[i]
	#	idx_tmp = length(idx_loc)==1 ? idx_loc[1] : idx_loc;
	#	price+=sum(weight_*mean(payoff(S[idx_tmp],abstractPayoff,rfCurve,maxT)) for (abstractPayoff,weight_) in dict_[i])
	#end

	price=sum( sum(weight_*mean(payoff(S[idx_[i]],abstractPayoff,rfCurve,maxT)) for (abstractPayoff,weight_) in dict_[i]) for i in IND_)
	
	return price;
end


function pricer(mcProcess::Dict{String,FinancialMonteCarlo.AbstractMonteCarloProcess},rfCurve::AbstractZeroRateCurve,mcConfig::MonteCarloConfiguration,dict_::Dict{String,Dict{FinancialMonteCarlo.AbstractPayoff,Number}})
	set_seed(mcConfig)
	underlyings_payoff=keys(dict_);
	underlyings_payoff_cpl=complete_2(underlyings_payoff,keys(mcProcess));
	price=0.0;
	for under_cpl in unique(underlyings_payoff_cpl)
		#options=dict_[under_]
		options=extract_(under_cpl,dict_)
		model=mcProcess[under_cpl]
		price=price+pricer(model,rfCurve,mcConfig,options);
	end
	
	#price= sum(pricer(mcProcess[under_],rfCurve,mcConfig,dict_[under_]) for under_ in underlyings_payoff);
	
	return price;
end

pricer_macro_array(BaseProcess)
pricer_macro_dict(BaseProcess)