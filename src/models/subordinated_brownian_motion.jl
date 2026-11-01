"""
    SubordinatedBrownianMotion(σ, drift, θ, subordinator)

Engine simulating a Brownian motion with drift `θ` and volatility `σ` evaluated at a subordinator,
plus a deterministic `drift`: at every time step
``X_{t+Δt} = X_t + \\text{drift}\\,Δt + θ\\,Δτ + σ \\sqrt{Δτ}\\,Z``, with ``Δτ`` drawn from `subordinator`.
If `drift` is a curve a [`FinancialMonteCarlo.SubordinatedBrownianMotionVec`](@ref) is returned.

# Arguments
- `σ`: volatility, must be positive.
- `drift`: deterministic drift in calendar time.
- `θ`: drift in subordinated time.
- `subordinator`: univariate continuous `Distribution` of the subordinator increment over one time step.
"""
struct SubordinatedBrownianMotion{num <: Number, num1 <: Number, num2 <: Number, Distr <: Distribution{Univariate, Continuous}, numtype <: Number} <: AbstractMonteCarloEngine{numtype}
    sigma::num
    drift::num1
    θ::num2
    subordinator_::Distr
    function SubordinatedBrownianMotion(sigma::num, drift::num1, θ::num2, dist::Distr) where {num <: Number, num1 <: Number, num2 <: Number, Distr <: Distribution{Univariate, Continuous}}
        ChainRulesCore.@ignore_derivatives @assert sigma > 0 "volatility must be positive"
        zero_typed = ChainRulesCore.@ignore_derivatives zero(num) + zero(num1) + zero(num2)
        return new{num, num1, num2, Distr, typeof(zero_typed)}(sigma, drift, θ, dist)
    end
end

export SubordinatedBrownianMotion;
