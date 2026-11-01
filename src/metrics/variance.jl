"""
    variance(mcProcess, rfCurve, mcConfig, payoff)

Sample variance of the discounted payoff over the simulated paths (a vector for a vector of payoffs or a smile).
The variance of the price estimate is this value divided by the number of paths.

# Arguments
- `mcProcess`: the model.
- `rfCurve`: the zero rate curve.
- `mcConfig`: the [`MonteCarloConfiguration`](@ref).
- `payoff`: a payoff, a vector of payoffs or a smile.
"""
function variance(mcProcess::BaseProcess, rfCurve::AbstractZeroRateCurve, mcConfig::MonteCarloConfiguration, abstractPayoff::AbstractPayoff)
    set_seed!(mcConfig)
    T = maturity(abstractPayoff)
    S = simulate(mcProcess, rfCurve, mcConfig, T)
    Payoff = payoff(S, abstractPayoff, rfCurve, mcConfig)
    variance_ = var(Payoff)
    return variance_
end

function variance(mcProcess::BaseProcess, rfCurve::AbstractZeroRateCurve, mcConfig::MonteCarloConfiguration, abstractPayoffs::Array{abstractPayoff_}) where {abstractPayoff_ <: AbstractPayoff}
    set_seed!(mcConfig)
    maxT = maximum([maturity(abstractPayoff) for abstractPayoff in abstractPayoffs])
    S = simulate(mcProcess, rfCurve, mcConfig, maxT)
    variance_ = [var(payoff(S, abstractPayoff, rfCurve, mcConfig, maxT)) for abstractPayoff in abstractPayoffs]

    return variance_
end

function variance(mcProcess::BaseProcess, rfCurve::AbstractZeroRateCurve, mcConfig::MonteCarloConfiguration, smile::AbstractSmile)
    set_seed!(mcConfig)
    T = maturity(smile)
    S = simulate(mcProcess, rfCurve, mcConfig, T)
    variance_ = vec(var(payoff(S, smile, rfCurve, mcConfig), dims = 1))

    return variance_
end
