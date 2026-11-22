"""
    delta(mcProcess, rfCurve, mcConfig, payoff, dS0 = 1e-7)
    delta(market, rfCurve, mcConfig, portfolio, underlying_name, dS0 = 1e-7)

Delta of `payoff` by forward finite difference: the spot of the underlying is bumped by `dS0` and
the payoff repriced with the same seed. For a [`MarketDataSet`](@ref) and a
[`FinancialMonteCarlo.Portfolio`](@ref), only the underlying called `underlying_name` is bumped
(also inside multivariate models); the delta is zero if it is not in the market data set.

# Arguments
- `mcProcess`: the model (a single name or multivariate process, or a [`MarketDataSet`](@ref)).
- `rfCurve`: the zero rate curve, e.g. [`ZeroRate`](@ref).
- `mcConfig`: the [`MonteCarloConfiguration`](@ref).
- `payoff`: a payoff, a vector of payoffs priced on the same paths, a smile, a [`FinancialMonteCarlo.Position`](@ref)
  or, with a [`MarketDataSet`](@ref), a [`FinancialMonteCarlo.Portfolio`](@ref).
- `dS0`: spot bump.

For exact pathwise deltas use a dual number as spot, see the Automatic differentiation section of the manual.
"""
function delta(mcProcess::BaseProcess, rfCurve::AbstractZeroRateCurve, mcConfig::MonteCarloConfiguration, abstractPayoffs::Union{AbstractPayoff, Array{<:AbstractPayoff}, AbstractSmile, Position}, dS0::Real = 1e-7)
    prices = pricer(mcProcess, rfCurve, mcConfig, abstractPayoffs)
    mcProcess_up = @set mcProcess.underlying.S0 += dS0
    prices_up = pricer(mcProcess_up, rfCurve, mcConfig, abstractPayoffs)
    return @. (prices_up - prices) / dS0
end

function delta(market::MarketDataSet, rfCurve::AbstractZeroRateCurve, mcConfig::MonteCarloConfiguration, portfolio::Portfolio, name::String, dS::Real = 1e-7)
    occursin('_', name) && throw(ArgumentError("deltas are defined on a single name, not $name"))
    price0 = pricer(market, rfCurve, mcConfig, portfolio)
    market_up = bump_spot(market, name, dS)
    isnothing(market_up) && return zero(price0)
    return (pricer(market_up, rfCurve, mcConfig, portfolio) - price0) / dS
end

# The market data set with the spot of name bumped by dS, alone or as a component of a multivariate model; nothing if
# name is not in it
function bump_spot(market::MarketDataSet, name::String, dS::Real)
    for (key, model) in market.models
        component = findfirst(==(name), underlying_names(key))
        isnothing(component) && continue
        models = copy(market.models)
        models[key] = bump_spot(model, component, dS)
        return MarketDataSet(models)
    end
    return nothing
end
bump_spot(model::AbstractMonteCarloProcess, ::Integer, dS::Real) = @set model.underlying.S0 += dS
function bump_spot(model::NDimensionalMonteCarloProcess, component::Integer, dS::Real)
    models = model.models
    bumped = @set models[component].underlying.S0 += dS
    return @set model.models = bumped
end
