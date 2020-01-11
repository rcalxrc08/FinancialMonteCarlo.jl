"""
General Interface for Pricing

		Delta=delta(mcProcess,rfCurve,mcBaseData,payoff_)
	
Where:\n
		mcProcess          = Process to be simulated.
		rfCurve  = Datas of the Spot.
		mcBaseData = Basic properties of MonteCarlo simulation
		payoff_ = Payoff(s) to be priced
		

		Delta     = Delta of the derivative

"""	
function delta(mcProcess::BaseProcess,rfCurve::AbstractZeroRateCurve,mcConfig::MonteCarloConfiguration,abstractPayoff::AbstractPayoff,dS0::Real=1e-7)

	Price=pricer(mcProcess,rfCurve,mcConfig,abstractPayoff);
	mcProcess_up=deepcopy(mcProcess)
	mcProcess_up.underlying.S0+=dS0;
	PriceUp=pricer(mcProcess_up,rfCurve,mcConfig,abstractPayoff);
	Delta=(PriceUp-Price)/dS0;

	return Delta;
end

function delta(mcProcess::BaseProcess,rfCurve::AbstractZeroRateCurve,mcConfig::MonteCarloConfiguration,abstractPayoffs::Array{AbstractPayoff},dS0::Real=1e-7)
		Prices=pricer(mcProcess,rfCurve,mcConfig,abstractPayoffs);
		mcProcess_up=deepcopy(mcProcess)
		mcProcess_up.underlying.S0+=dS0;
		PricesUp=pricer(mcProcess_up,spotData_1,mcConfig,abstractPayoffs);
		Delta=(PricesUp.-Prices)./dS0;
	
	return Delta;
end

function delta(mcProcess::BaseProcess,rfCurve::AbstractZeroRateCurve,mcConfig::MonteCarloConfiguration,abstractPayoffs::Dict{FinancialMonteCarlo.AbstractPayoff,Number},dS0::Real=1e-7)
		Prices=pricer(mcProcess,rfCurve,mcConfig,abstractPayoffs);
		mcProcess_up=deepcopy(mcProcess)
		mcProcess_up.underlying.S0+=dS0;
		PricesUp=pricer(mcProcess_up,rfCurve,mcConfig,abstractPayoffs);
		Delta=(PricesUp.-Prices)./dS0;
	
	return Delta;
end

function delta(mcProcess::Dict{String,FinancialMonteCarlo.AbstractMonteCarloProcess},rfCurve::AbstractZeroRateCurve,mcConfig::MonteCarloConfiguration,dict_::Dict{String,Dict{FinancialMonteCarlo.AbstractPayoff,Number}},x::String)
	set_seed(mcConfig)
	underlyings_models=keys(mcProcess)
	underlyings_payoff=keys(dict_)
	price0=pricer(mcProcess,rfCurve,mcConfig,dict_);
	delta_=0.0;
	deps_=derive_dep(x,mcProcess);
	mcProcess_=mcProcess;
	for under_ in underlyings_payoff
		options=dict_[under_]
		model=mcProcess[under_]
		price=price+pricer(model,rfCurve,mcConfig,options);
	end
	
	return price;
end