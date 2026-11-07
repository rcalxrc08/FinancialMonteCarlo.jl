include("utils.jl")

"""
    pricer(mcProcess, rfCurve, mcConfig, payoff)

Monte Carlo estimate of the price of `payoff`: the mean of the discounted payoff over the
simulated paths. The random number generator is reset to the seed of `mcConfig` first, so the
result is reproducible. Returns a vector of prices for a vector of payoffs or a smile; for a smile
in `CudaMode` it is a `CuArray`, so the computations on the prices can stay on the GPU.

# Arguments
- `mcProcess`: the model (a single name or multivariate process, or a [`MarketDataSet`](@ref)).
- `rfCurve`: the zero rate curve, e.g. [`ZeroRate`](@ref).
- `mcConfig`: the [`MonteCarloConfiguration`](@ref).
- `payoff`: a payoff, a vector of payoffs priced on the same paths, a smile, a [`FinancialMonteCarlo.Position`](@ref)
  or, with a [`MarketDataSet`](@ref), a [`FinancialMonteCarlo.Portfolio`](@ref).

# Examples
```julia
model = BlackScholesProcess(0.2, Underlying(100.0, 0.01))
mcConfig = MonteCarloConfiguration(10_000, 30)
pricer(model, ZeroRate(0.02), mcConfig, EuropeanOption(1.0, 100.0))
pricer(model, ZeroRate(0.02), mcConfig, [Forward(1.0), AmericanOption(1.0, 100.0, false)])
```
"""
function pricer(mcProcess::BaseProcess, rfCurve::AbstractZeroRateCurve, mcConfig::AbstractMonteCarloConfiguration, abstractPayoff::AbstractPayoff)
    set_seed!(mcConfig)
    T = maturity(abstractPayoff)
    S = simulate(mcProcess, rfCurve, mcConfig, T)
    Payoff = payoff(S, abstractPayoff, rfCurve, mcConfig)
    Price = mean(Payoff)
    return Price
end

function pricer(mcProcess::BaseProcess, rfCurve::AbstractZeroRateCurve, mcConfig::AbstractMonteCarloConfiguration, smile::AbstractSmile)
    abstractPayoff = to_mode(smile, execution_mode(mcConfig))
    set_seed!(mcConfig)
    T = maturity(abstractPayoff)
    S = simulate(mcProcess, rfCurve, mcConfig, T)
    return convert(AbstractArray{price_type(mcProcess, rfCurve, abstractPayoff)}, vec(mean(payoff(S, abstractPayoff, rfCurve, mcConfig, T), dims = 1)))
end

# The arrays of the paths with numbers of type T (see price_type), on the device of the mode
get_matrix_type(mcConfig::MonteCarloConfiguration, ::BaseProcess, ::Type{T}) where {T <: Number} = array_type(mcConfig.mode){T}(undef, mcConfig.Nsim, mcConfig.Nstep + 1);
get_array_type(mcConfig::MonteCarloConfiguration, ::BaseProcess, ::Type{T}) where {T <: Number} = array_type(mcConfig.mode){T}(undef, mcConfig.Nsim);
get_matrix_type(mcConfig::MonteCarloConfiguration, ::VectorialMonteCarloProcess, ::Type{T}) where {T <: Number} = Array{array_type(mcConfig.mode){T, 2}};

function pricer(mcProcess::BaseProcess, rfCurve::AbstractZeroRateCurve, mcConfig::MonteCarloConfiguration, abstractPayoffs::Array{abstractPayoff_, N}) where {abstractPayoff_ <: AbstractPayoff, N}
    set_seed!(mcConfig)
    maxT = maturity(abstractPayoffs)
    S = simulate(mcProcess, rfCurve, mcConfig, maxT)
    T_price = price_type(mcProcess, rfCurve, abstractPayoffs)
    Prices::Array{T_price, N} = [mean(payoff(S, abstractPayoff, rfCurve, mcConfig, maxT)) for abstractPayoff in abstractPayoffs]

    return Prices
end

function pricer(mcProcess::BaseProcess, rfCurve::AbstractZeroRateCurve, mcConfig::MonteCarloConfiguration, dict_::Position)
    set_seed!(mcConfig)
    maxT = maturity(dict_)
    S = simulate(mcProcess, rfCurve, mcConfig, maxT)
    price = sum(weight_ * mean(payoff(S, abstractPayoff, rfCurve, mcConfig, maxT)) for (abstractPayoff, weight_) in dict_)

    return price
end

#####Pricer for multivariate
function pricer(mcProcess::VectorialMonteCarloProcess, rfCurve::AbstractZeroRateCurve, mcConfig::MonteCarloConfiguration, special_array_payoff::Array{Position})
    set_seed!(mcConfig)
    N_ = length(mcProcess.models)
    idx_ = compute_indices(N_)
    ## Filter special_array_payoff for undef
    IND_ = collect(1:length(special_array_payoff))
    filter!(i -> isassigned(special_array_payoff, i), IND_)
    dict_cl = special_array_payoff[IND_]
    #Compute Maximum time to maturity, all of the underlyings will be simulated till maxT.
    maxT = maturity(dict_cl)
    S = simulate(mcProcess, rfCurve, mcConfig, maxT)
    price = sum(sum(weight_ * mean(payoff(S[idx_[i]], abstractPayoff, rfCurve, mcConfig, maxT)) for (abstractPayoff, weight_) in special_array_payoff[i]) for i in IND_)

    return price
end

function pricer(mcProcess::MarketDataSet, rfCurve::AbstractZeroRateCurve, mcConfig::MonteCarloConfiguration, dict_::Portfolio)
    set_seed!(mcConfig)
    underlyings_payoff = keys(dict_)
    underlyings_payoff_cpl = get_underlyings_identifier(underlyings_payoff, keys(mcProcess))
    T_price = price_type(rfCurve, collect(values(dict_)), collect(values(mcProcess)))
    price::T_price = sum(pricer(mcProcess[under_cpl], rfCurve, mcConfig, extract_option_from_portfolio(under_cpl, dict_)) for under_cpl in unique(underlyings_payoff_cpl))

    return price
end
