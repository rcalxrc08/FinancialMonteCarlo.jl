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
    return mean(payoff(S, abstractPayoff, rfCurve, mcConfig))
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
    prices::Array{T_price, N} = [mean(payoff(S, abstractPayoff, rfCurve, mcConfig, maxT)) for abstractPayoff in abstractPayoffs]

    return prices
end

# A position is priced path by path: the mean of its discounted payoff, the weighted sum of those of its legs on the same
# paths (so that its variance and confidence interval are those of the position, see payoff_moments)
function pricer(mcProcess::BaseProcess, rfCurve::AbstractZeroRateCurve, mcConfig::MonteCarloConfiguration, position::Position)
    T_price = price_type(mcProcess, rfCurve, position)
    isempty(position) && return zero(T_price)
    set_seed!(mcConfig)
    maxT = maturity(position)
    S = simulate(mcProcess, rfCurve, mcConfig, maxT)
    price::T_price = mean(payoff(S, position, rfCurve, mcConfig, maxT))
    return price
end

#####Pricer for multivariate
# One position per subset of the components of the model, in the order of compute_indices; undefined or empty where
# there is none. The paths of a subset are those of its components.
function pricer(mcProcess::VectorialMonteCarloProcess, rfCurve::AbstractZeroRateCurve, mcConfig::MonteCarloConfiguration, positions::Array{<:Position})
    subsets = compute_indices(length(mcProcess.models))
    priced = [i for i in eachindex(positions) if isassigned(positions, i) && !isempty(positions[i])]
    T_price = price_type(mcProcess, rfCurve, positions[priced])
    isempty(priced) && return zero(T_price)
    set_seed!(mcConfig)
    maxT = maximum(i -> maturity(positions[i]), priced)
    S = simulate(mcProcess, rfCurve, mcConfig, maxT)
    price::T_price = sum(mean(payoff(S[subsets[i]], positions[i], rfCurve, mcConfig, maxT)) for i in priced)
    return price
end

# Each model of the market data set prices the positions on its underlyings (for a multivariate model, also those on its
# components and their subsets). Positions on names missing from the market data set are not priced.
# The generator is seeded once: the models are simulated one after the other, in the order of the market data set, each
# continuing the random numbers of the previous one, so that the underlyings are independent (#38). Reseeding for each
# model would give them all the same numbers, i.e. perfectly correlated paths.
function pricer(market::MarketDataSet, rfCurve::AbstractZeroRateCurve, mcConfig::MonteCarloConfiguration, portfolio::Portfolio)
    T_price = price_type(rfCurve, portfolio, market)
    priced = [(model, positions) for (model, positions) in ((model, positions_on(portfolio, key)) for (key, model) in market.models) if !isnothing(positions)]
    check_independent_streams(mcConfig.sampling, length(priced))
    set_seed!(mcConfig)
    continuing = @set mcConfig.reseed = false
    price::T_price = zero(T_price)
    for (model, positions) in priced
        price += pricer(model, rfCurve, continuing, positions)
    end
    return price
end

# A Sobol sequence does not depend on the generator: every model would get the same points
check_independent_streams(::AbstractSampling, ::Integer) = nothing
check_independent_streams(::SobolMC, n::Integer) = ChainRulesCore.@ignore_derivatives n <= 1 || throw(ArgumentError("SobolMC gives every model of a market data set the same points: price a portfolio on several models with another sampling"))
