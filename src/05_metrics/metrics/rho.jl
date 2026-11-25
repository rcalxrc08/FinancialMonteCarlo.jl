bump_rate(rfCurve::ZeroRate{<:Number}, dr) = ZeroRate(rfCurve.r + dr)
bump_rate(rfCurve::ZeroRate{<:Curve}, dr) = ZeroRate(map(x -> x + dr, rfCurve.r))

"""
    rho(mcProcess, underlying, rfCurve, mcConfig, payoff, dr = 1e-7)

Rho of `payoff` by forward finite difference: the zero rate (every node of a
curve of a [`ZeroRate`](@ref), i.e. a parallel shift) is bumped by `dr` and the payoff repriced with the
same seed. Returns a vector for a vector of payoffs or a smile.

# Arguments
- `mcProcess`: the model.
- `underlying`: the [`Underlying`](@ref) (spot and dividends), a tuple of them for a multivariate model.
- `rfCurve`: the zero rate curve.
- `mcConfig`: the [`MonteCarloConfiguration`](@ref).
- `payoff`: a payoff, a vector of payoffs, a smile or a [`FinancialMonteCarlo.Position`](@ref).
- `dr`: rate bump.
"""
function rho(mcProcess::BaseProcess, und::UnderlyingArg, rfCurve::AbstractZeroRateCurve, mcConfig::MonteCarloConfiguration, abstractPayoffs::Union{AbstractPayoff, Array{<:AbstractPayoff}, AbstractSmile, Position}, dr::Real = 1e-7)
    prices = pricer(mcProcess, und, rfCurve, mcConfig, abstractPayoffs)
    prices_up = pricer(mcProcess, und, bump_rate(rfCurve, dr), mcConfig, abstractPayoffs)
    return @. (prices_up - prices) / dr
end

export rho;
