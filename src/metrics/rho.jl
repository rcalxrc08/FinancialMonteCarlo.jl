bump_rate(rfCurve::ZeroRate, dr) = ZeroRate(rfCurve.r + dr)
bump_rate(rfCurve::ZeroRateCurve, dr) = ZeroRateCurve(map(x -> x + dr, rfCurve.r))

"""
    rho(mcProcess, rfCurve, mcConfig, payoff, dr = 1e-7)

Rho of `payoff` by forward finite difference: the zero rate (every node of a
[`ZeroRateCurve`](@ref), i.e. a parallel shift) is bumped by `dr` and the payoff repriced with the
same seed. Returns a vector for a vector of payoffs or a smile.

# Arguments
- `mcProcess`: the model.
- `rfCurve`: the zero rate curve.
- `mcConfig`: the [`MonteCarloConfiguration`](@ref).
- `payoff`: a payoff, a vector of payoffs or a smile.
- `dr`: rate bump.
"""
function rho(mcProcess::BaseProcess, rfCurve::AbstractZeroRateCurve, mcConfig::MonteCarloConfiguration, abstractPayoffs::Union{AbstractPayoff, Array{<:AbstractPayoff}, AbstractSmile}, dr::Real = 1e-7)
    Prices = pricer(mcProcess, rfCurve, mcConfig, abstractPayoffs)
    PricesUp = pricer(mcProcess, bump_rate(rfCurve, dr), mcConfig, abstractPayoffs)
    rho = @. (PricesUp - Prices) / dr

    return rho
end

export rho;
