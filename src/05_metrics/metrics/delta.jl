"""
    delta(mcProcess, underlying, rfCurve, mcConfig, payoff, dS0 = 1e-7)
    delta(market, rfCurve, mcConfig, portfolio, underlying_name, dS0 = 1e-7)

Delta of `payoff` by forward finite difference: the spot of the underlying is bumped by `dS0` and
the payoff repriced with the same seed. For a [`MarketDataSet`](@ref) and a
[`FinancialMonteCarlo.Portfolio`](@ref), only the underlying called `underlying_name` is bumped
(also inside multivariate models); the delta is zero if it is not in the market data set.

# Arguments
- `mcProcess`: a single name model, or a [`MarketDataSet`](@ref) without `underlying`.
- `underlying`: the [`Underlying`](@ref) whose spot is bumped.
- `rfCurve`: the zero rate curve, e.g. [`ZeroRate`](@ref).
- `mcConfig`: the [`MonteCarloConfiguration`](@ref).
- `payoff`: a payoff, a vector of payoffs priced on the same paths, a smile, a [`FinancialMonteCarlo.Position`](@ref)
  or, with a [`MarketDataSet`](@ref), a [`FinancialMonteCarlo.Portfolio`](@ref).
- `dS0`: spot bump.

For exact pathwise deltas use a dual number as spot, see the Automatic differentiation section of the manual.
"""
function delta(mcProcess::BaseProcess, und::AbstractUnderlying, rfCurve::AbstractZeroRateCurve, mcConfig::MonteCarloConfiguration, abstractPayoffs::Union{AbstractPayoff, Array{<:AbstractPayoff}, AbstractSmile, Position}, dS0::Real = 1e-7)
    prices = pricer(mcProcess, und, rfCurve, mcConfig, abstractPayoffs)
    und_up = bump_spot(und, dS0)
    prices_up = pricer(mcProcess, und_up, rfCurve, mcConfig, abstractPayoffs)
    return @. (prices_up - prices) / dS0
end

function delta(market::MarketDataSet, rfCurve::AbstractZeroRateCurve, mcConfig::MonteCarloConfiguration, portfolio::Portfolio, name::String, dS::Real = 1e-7)
    occursin('_', name) && throw(ArgumentError("deltas are defined on a single name, not $name"))
    price0 = pricer(market, rfCurve, mcConfig, portfolio)
    market_up = bump_spot(market, name, dS)
    isnothing(market_up) && return zero(price0)
    return (pricer(market_up, rfCurve, mcConfig, portfolio) - price0) / dS
end

# The market data set with the spot of name bumped by dS: its underlying, shared by every model of the name (alone or
# as a component of a multivariate model); nothing if name has no underlying in it
function bump_spot(market::MarketDataSet, name::String, dS::Real)
    haskey(market.underlyings, name) || return nothing
    underlyings = copy(market.underlyings)
    underlyings[name] = bump_spot(underlyings[name], dS)
    return MarketDataSet(underlyings, market.models)
end
bump_spot(und::AbstractUnderlying, dS::Real) = @set und.S0 += dS
