"""
    BrownianMotionVec(σ, μ)

Engine simulating a Brownian motion with time dependent drift, ``X_t = \\int_0^t μ(s)\\,ds + σ W_t``.
Built by [`BrownianMotion`](@ref) when the drift is a curve.

# Arguments
- `σ`: volatility, must be positive.
- `μ`: drift curve.
"""
struct BrownianMotionVec{num <: Number, num1 <: Number, num4 <: Number, numtype <: Number} <: AbstractMonteCarloEngine{numtype}
    σ::num
    μ::CurveType{num1, num4}
    function BrownianMotionVec(σ::num, μ::CurveType{num1, num4}) where {num <: Number, num1 <: Number, num4 <: Number}
        ChainRulesCore.@ignore_derivatives check_domain(σ > 0, σ, "Volatility must be positive")
        zero_typed = ChainRulesCore.@ignore_derivatives zero(num) + zero(num1) + zero(num4)
        return new{num, num1, num4, typeof(zero_typed)}(σ, μ)
    end
end

function BrownianMotion(σ::num, μ::CurveType{num1, num4}) where {num <: Number, num1 <: Number, num4 <: Number}
    return BrownianMotionVec(σ, μ)
end

function simulate!(X, mcProcess::Union{BrownianMotion, BrownianMotionVec}, mcBaseData::MonteCarloConfiguration, T::Number)
    σ = mcProcess.σ
    μ = mcProcess.μ
    ChainRulesCore.@ignore_derivatives check_domain(T > 0, T, "Time to Maturity must be positive")
    Nstep = mcBaseData.Nstep
    dt = T / Nstep
    stddev_bm = σ * sqrt(dt)
    drifts = step_drifts(μ, dt, Nstep)
    isDualZero = stddev_bm * 0 * drift_of_step(drifts, 1)
    fill!(view(X, :, 1), isDualZero)
    sides = path_sides(X, mcBaseData)
    Z = draw_buffer(mcBaseData, typeof(get_rng_type(isDualZero)), draws_per_side(sides))
    @inbounds for j = 1:Nstep
        drift = drift_of_step(drifts, j)
        draw_normals!(Z, mcBaseData)
        foreach(sides) do (X_, sign_)
            bm_step!(view(X_, :, j + 1), view(X_, :, j), drift, sign_ * stddev_bm, Z)
        end
    end
    return
end
