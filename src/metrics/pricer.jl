include("utils.jl")

"""
Low Level Interface for Pricing

		Price=pricer(mcProcess,rfCurve,mcConfig,abstractPayoff)
	
Where:\n
		mcProcess          = Process to be simulated.
		rfCurve  = Zero Rate Data.
		mcConfig = Basic properties of MonteCarlo simulation
		abstractPayoff = Payoff(s) to be priced

		Price     = Price of the derivative

"""
function pricer(mcProcess::BaseProcess, rfCurve::AbstractZeroRateCurve, mcConfig::AbstractMonteCarloConfiguration, abstractPayoff::AbstractPayoff)
    set_seed(mcConfig)
    T = maturity(abstractPayoff)
    S = simulate(mcProcess, rfCurve, mcConfig, T)
    Payoff = payoff(S, abstractPayoff, rfCurve, mcConfig)
    Price = mean(Payoff)
    return Price
end

get_matrix_type(mcConfig::MonteCarloConfiguration{<:Integer, <:Integer, <:AbstractMonteCarloMethod, <:BaseMode, <:Random.AbstractRNG}, ::BaseProcess, price) = Matrix{typeof(price)}(undef, mcConfig.Nsim, mcConfig.Nstep + 1);
get_array_type(mcConfig::MonteCarloConfiguration{<:Integer, <:Integer, <:AbstractMonteCarloMethod, <:BaseMode, <:Random.AbstractRNG}, ::BaseProcess, price) = Array{typeof(price)}(undef, mcConfig.Nstep);
get_matrix_type(::MonteCarloConfiguration{<:Integer, <:Integer, <:AbstractMonteCarloMethod, <:BaseMode, <:Random.AbstractRNG}, ::VectorialMonteCarloProcess, price) = Array{Matrix{typeof(price)}};

function pricer(mcProcess::BaseProcess, rfCurve::AbstractZeroRateCurve, mcConfig::MonteCarloConfiguration, abstractPayoffs::Array{abstractPayoff_}) where {abstractPayoff_ <: AbstractPayoff}
    set_seed(mcConfig)
    maxT = maximum(maturity.(abstractPayoffs))
    S = simulate(mcProcess, rfCurve, mcConfig, maxT)
    zero_typed = predict_output_type_zero(mcProcess, rfCurve, mcConfig, abstractPayoffs)
    Prices::Array{typeof(zero_typed)} = [mean(payoff(S, abstractPayoff, rfCurve, mcConfig, maxT)) for abstractPayoff in abstractPayoffs]

    return Prices
end

function pricer(mcProcess::BaseProcess, rfCurve::AbstractZeroRateCurve, mcConfig::MonteCarloConfiguration, dict_::Position)
    set_seed(mcConfig)
    abstractPayoffs = keys(dict_)
    maxT = maximum([maturity(abstractPayoff) for abstractPayoff in abstractPayoffs])
    S = simulate(mcProcess, rfCurve, mcConfig, maxT)
    price = sum(weight_ * mean(payoff(S, abstractPayoff, rfCurve, mcConfig, maxT)) for (abstractPayoff, weight_) in dict_)

    return price
end

#####Pricer for multivariate
function pricer(mcProcess::VectorialMonteCarloProcess, rfCurve::AbstractZeroRateCurve, mcConfig::MonteCarloConfiguration, special_array_payoff::Array{Position})
    set_seed(mcConfig)
    N_ = length(mcProcess.models)
    idx_ = compute_indices(N_)
    ## Filter special_array_payoff for undef
    IND_ = collect(1:length(special_array_payoff))
    filter!(i -> isassigned(special_array_payoff, i), IND_)
    dict_cl = special_array_payoff[IND_]
    #Compute Maximum time to maturity, all of the underlyings will be simulated till maxT.
    maxT = maximum([maximum(maturity.(collect(keys(ar_el)))) for ar_el in dict_cl])
    S = simulate(mcProcess, rfCurve, mcConfig, maxT)
    #price=0.0;
    #for i in IND_
    #	idx_loc=idx_[i]
    #	idx_tmp = length(idx_loc)==1 ? idx_loc[1] : idx_loc;
    #	price+=sum(weight_*mean(payoff(S[idx_tmp],abstractPayoff,rfCurve,maxT)) for (abstractPayoff,weight_) in special_array_payoff[i])
    #end

    price = sum(sum(weight_ * mean(payoff(S[idx_[i]], abstractPayoff, rfCurve, mcConfig, maxT)) for (abstractPayoff, weight_) in special_array_payoff[i]) for i in IND_)

    return price
end

function pricer(mcProcess::MarketDataSet, rfCurve::AbstractZeroRateCurve, mcConfig::MonteCarloConfiguration, dict_::Portfolio)
    set_seed(mcConfig)
    underlyings_payoff = keys(dict_)
    underlyings_payoff_cpl = get_underlyings_identifier(underlyings_payoff, keys(mcProcess))
    #price=0.0;
    #for under_cpl in unique(underlyings_payoff_cpl)
    #	#options=dict_[under_]
    #	options=extract_option_from_portfolio(under_cpl,dict_)
    #	model=mcProcess[under_cpl]
    #	price=price+pricer(model,rfCurve,mcConfig,options);
    #end
    zero_typed = predict_output_type_zero(rfCurve, mcConfig, collect(keys(collect(values(dict_)))), collect(values(mcProcess)))
    price::typeof(zero_typed) = sum(pricer(mcProcess[under_cpl], rfCurve, mcConfig, extract_option_from_portfolio(under_cpl, dict_)) for under_cpl in unique(underlyings_payoff_cpl))

    return price
end
