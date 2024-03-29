"""
General Interface for Delta

		Delta=delta(mcProcess,rfCurve,mcConfig,abstractPayoff,dS0)
	
Where:\n
		mcProcess      = Process to be simulated.
		rfCurve        = Zero Rate Data.
		mcConfig       = Basic properties of MonteCarlo simulation
		abstractPayoff = Payoff(s) to be priced
		dS0            = increment [optional, default to 1e-7]
		
		Delta          = Delta of the derivative

"""
function delta(mcProcess::BaseProcess, rfCurve::AbstractZeroRateCurve, mcConfig::MonteCarloConfiguration, abstractPayoff::AbstractPayoff, dS0::Real = 1e-7)
    Price = pricer(mcProcess, rfCurve, mcConfig, abstractPayoff)
    mcProcess_up = @set mcProcess.underlying.S0 += dS0
    set_seed!(mcConfig)
    PriceUp = pricer(mcProcess_up, rfCurve, mcConfig, abstractPayoff)
    Delta = (PriceUp - Price) / dS0

    return Delta
end

function delta(mcProcess::BaseProcess, rfCurve::AbstractZeroRateCurve, mcConfig::MonteCarloConfiguration, abstractPayoffs::Array{abstractPayoff}, dS0::Real = 1e-7) where {abstractPayoff <: AbstractPayoff}
    Prices = pricer(mcProcess, rfCurve, mcConfig, abstractPayoffs)
    mcProcess_up = @set mcProcess.underlying.S0 += dS0
    PricesUp = pricer(mcProcess_up, rfCurve, mcConfig, abstractPayoffs)
    Delta = @. (PricesUp - Prices) / dS0

    return Delta
end

function delta(mcProcess::BaseProcess, rfCurve::AbstractZeroRateCurve, mcConfig::MonteCarloConfiguration, abstractPayoffs::Dict{AbstractPayoff, Number}, dS0::Real = 1e-7)
    Prices = pricer(mcProcess, rfCurve, mcConfig, abstractPayoffs)
    mcProcess_up = @set mcProcess.underlying.S0 += dS0
    set_seed!(mcConfig)
    PricesUp = pricer(mcProcess_up, rfCurve, mcConfig, abstractPayoffs)
    Delta = @. (PricesUp - Prices) / dS0

    return Delta
end

function delta(mcProcess::Dict{String, AbstractMonteCarloProcess}, rfCurve::AbstractZeroRateCurve, mcConfig::MonteCarloConfiguration, dict_::Dict{String, Dict{AbstractPayoff, Number}}, underl_::String, dS::Real = 1e-7)
    ChainRulesCore.@ignore_derivatives @assert isnothing(findfirst("_", underl_)) "deltas are defined on single name"
    set_seed!(mcConfig)
    price0 = pricer(mcProcess, rfCurve, mcConfig, dict_)
    price2 = 0.0
    set_seed!(mcConfig)
    mcProcess_up = deepcopy(mcProcess)
    keys_mkt = collect(keys(mcProcess_up))
    idx_supp = 0
    idx_1 = findfirst(y -> y == underl_, keys_mkt)
    if (isnothing(idx_1))
        idx_1 = findfirst(out_el -> any(out_el_sub -> out_el_sub == underl_, split(out_el, "_")), keys_mkt)
        if (isnothing(idx_1))
            return 0.0
        end
        multi_name = keys_mkt[idx_1]
        idx_supp = findfirst(x_ -> x_ == underl_, split(multi_name, "_"))
        model_ = mcProcess_up[multi_name]
        ok_m = model_.models[idx_supp]
        ok_m = @set ok_m.underlying.S0 += dS
        # model_.models[idx_supp].underlying.S0 += dS
        models_ = model_.models
        new_models = @set models_[idx_supp] = ok_m
        new_model = @set model_.models = new_models
        delete!(mcProcess_up, multi_name)
        tmp_mkt = multi_name → new_model
        mcProcess_up = mcProcess_up + tmp_mkt
        price2 = pricer(mcProcess_up, rfCurve, mcConfig, dict_)
    else
        model = mcProcess_up[keys_mkt[idx_1]]
        model_up = @set model.underlying.S0 += dS
        delete!(mcProcess_up, keys_mkt[idx_1])
        tmp_mkt = keys_mkt[idx_1] → model_up
        mcProcess_up = mcProcess_up + tmp_mkt
        price2 = pricer(mcProcess_up, rfCurve, mcConfig, dict_)
    end

    return (price2 - price0) / dS
end
