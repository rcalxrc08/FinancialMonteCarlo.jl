"""
    SubordinatedBrownianMotionVec(σ, drift, θ, subordinator)

[`SubordinatedBrownianMotion`](@ref) with a time dependent deterministic drift.

# Arguments
- `σ`: volatility, must be positive.
- `drift`: drift curve in calendar time.
- `θ`: drift in subordinated time.
- `subordinator`: univariate continuous `Distribution` of the subordinator increment over one time step.
"""
struct SubordinatedBrownianMotionVec{num <: Number, num1 <: Number, num4 <: Number, num2 <: Number, Distr <: Distribution{Univariate, Continuous}, numtype <: Number} <: AbstractMonteCarloEngine{numtype}
    sigma::num
    drift::CurveType{num1, num4}
    θ::num2
    subordinator_::Distr
    function SubordinatedBrownianMotionVec(sigma::num, drift::CurveType{num1, num4}, θ::num2, dist::Distr) where {num <: Number, num1 <: Number, num2 <: Number, num4 <: Number, Distr <: Distribution{Univariate, Continuous}}
        ChainRulesCore.@ignore_derivatives @assert sigma > 0 "volatility must be positive"
        zero_typed = ChainRulesCore.@ignore_derivatives zero(num) + zero(num1) + zero(num4) + zero(num2)
        return new{num, num1, num4, num2, Distr, typeof(zero_typed)}(sigma, drift, θ, dist)
    end
end

function SubordinatedBrownianMotion(σ::num, drift::CurveType{num1, num4}, θ::num2, subordinator_::Distr) where {num <: Number, num1 <: Number, num2 <: Number, num4 <: Number, Distr <: Distribution{Univariate, Continuous}}
    return SubordinatedBrownianMotionVec(σ, drift, θ, subordinator_)
end

const AnySubordinatedBrownianMotion = Union{SubordinatedBrownianMotion, SubordinatedBrownianMotionVec}

# The drifts, the type of the subordinator increments and the zero carrying the derivatives of the inputs
function subordinated_setup(mcProcess::AnySubordinatedBrownianMotion, mcBaseData::MonteCarloConfiguration, T::Number)
    ChainRulesCore.@ignore_derivatives @assert T > 0
    drifts = step_drifts(mcProcess.drift, T / mcBaseData.Nstep, mcBaseData.Nstep)
    type_sub = typeof(quantile(mcProcess.subordinator_, 1 // 2))
    return drifts, type_sub, mcProcess.sigma * zero(type_sub) * 0 * drift_of_step(drifts, 1)
end

function simulate!(X, mcProcess::AnySubordinatedBrownianMotion, mcBaseData::MonteCarloConfiguration, T::Number)
    sigma = mcProcess.sigma
    θ = mcProcess.θ
    drifts, type_sub, isDualZero = subordinated_setup(mcProcess, mcBaseData, T)
    fill!(view(X, :, 1), isDualZero)
    sides = path_sides(X, mcBaseData)
    N = draws_per_side(sides)
    dt_s = draw_buffer(mcBaseData, type_sub, N)
    dt_s_host = subordinator_host(dt_s)
    Z = draw_buffer(mcBaseData, typeof(get_rng_type(isDualZero)), N)
    @inbounds for j = 1:mcBaseData.Nstep
        draw_subordinator!(dt_s, dt_s_host, mcProcess.subordinator_, mcBaseData)
        draw_normals!(Z, mcBaseData)
        drift = drift_of_step(drifts, j)
        for (X_, sign_) in sides
            sigma_step = sign_ * sigma
            # SUBORDINATED BROWNIAN MOTION (dt_s=time change)
            @views @. X_[:, j+1] = X_[:, j] + drift + θ * dt_s + sigma_step * sqrt(dt_s) * Z
        end
    end
    return
end
