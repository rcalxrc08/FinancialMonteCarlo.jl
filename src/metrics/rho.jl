"""
    rho(mcProcess, rfCurve, mcConfig, payoff, dr = 1e-7)

Rho of `payoff` by forward finite difference: the zero rate (every node of a
[`ZeroRateCurve`](@ref), i.e. a parallel shift) is bumped by `dr` and the payoff repriced with the
same seed. Returns a vector for a vector of payoffs.

# Arguments
- `mcProcess`: the model.
- `rfCurve`: the zero rate curve.
- `mcConfig`: the [`MonteCarloConfiguration`](@ref).
- `payoff`: a payoff or a vector of payoffs.
- `dr`: rate bump.
"""
bump_rate(rfCurve::ZeroRate, dr) = ZeroRate(rfCurve.r + dr)
bump_rate(rfCurve::ZeroRateCurve, dr) = ZeroRateCurve(map(x -> x + dr, rfCurve.r))

function rho(mcProcess::BaseProcess, rfCurve::AbstractZeroRateCurve, mcConfig::MonteCarloConfiguration, abstractPayoff::AbstractPayoff, dr::Real = 1e-7)
    Price = pricer(mcProcess, rfCurve, mcConfig, abstractPayoff)
    rfCurve_1 = bump_rate(rfCurve, dr)
    PriceUp = pricer(mcProcess, rfCurve_1, mcConfig, abstractPayoff)
    rho = (PriceUp - Price) / dr

    return rho
end

function rho(mcProcess::BaseProcess, rfCurve::AbstractZeroRateCurve, mcConfig::MonteCarloConfiguration, abstractPayoffs::Array{abstractPayoff_}, dr::Real = 1e-7) where {abstractPayoff_ <: AbstractPayoff}
    Prices = pricer(mcProcess, rfCurve, mcConfig, abstractPayoffs)
    rfCurve_1 = bump_rate(rfCurve, dr)
    PricesUp = pricer(mcProcess, rfCurve_1, mcConfig, abstractPayoffs)
    rho = @. (PricesUp - Prices) / dr

    return rho
end

export rho;
