"""
Struct for Brownian Motion

		bmProcess=BrownianMotionVec(σ::num1,μ::num2) where {num1,num2 <: Number}
	
Where:

		σ	=	volatility of the process.
		μ	=	drift of the process.
"""
mutable struct BrownianMotionVec{num <: Number, num1 <: Number, num4 <: Number, numtype <: Number} <: AbstractMonteCarloEngine{numtype}
    σ::num
    μ::CurveType{num1, num4}
    function BrownianMotionVec(σ::num, μ::CurveType{num1, num4}) where {num <: Number, num1 <: Number, num4 <: Number}
        @assert σ > 0 "Volatility must be positive"
        zero_typed = zero(num) + zero(num1) + zero(num4)
        return new{num, num1, num4, typeof(zero_typed)}(σ, μ)
    end
end

function BrownianMotion(σ::num, μ::CurveType{num1, num4}) where {num <: Number, num1 <: Number, num4 <: Number}
    return BrownianMotionVec(σ, μ)
end

function simulate!(X, mcProcess::BrownianMotionVec, mcBaseData::SerialMonteCarloConfig, T::numb) where {numb <: Number}
    Nsim = mcBaseData.Nsim
    Nstep = mcBaseData.Nstep
    σ = mcProcess.σ
    μ = mcProcess.μ
    @assert T > 0
    dt = T / Nstep
    stddev_bm = σ * sqrt(dt)
    zero_drift = incremental_integral(μ, dt * 0, dt)
    isDualZero = stddev_bm * 0 * zero_drift
    view(X, :, 1) .= isDualZero
    Z = Array{typeof(get_rng_type(isDualZero))}(undef, Nsim)
    @inbounds for j = 1:Nstep
        tmp_ = incremental_integral(μ, (j - 1) * dt, dt)
        randn!(mcBaseData.parallelMode.rng, Z)
        @views @. X[:, j+1] = X[:, j] + tmp_ + stddev_bm * Z
    end

    nothing
end

function simulate!(X, mcProcess::BrownianMotionVec, mcBaseData::SerialAntitheticMonteCarloConfig, T::numb) where {numb <: Number}
    Nsim = mcBaseData.Nsim
    Nstep = mcBaseData.Nstep
    σ = mcProcess.σ
    μ = mcProcess.μ
    @assert T > 0
    dt = T / Nstep
    stddev_bm = σ * sqrt(dt)
    zero_drift = incremental_integral(μ, dt * 0, dt)
    isDualZero = stddev_bm * 0 * zero_drift
    view(X, :, 1) .= isDualZero
    Nsim_2 = div(Nsim, 2)
    Z = Array{typeof(get_rng_type(isDualZero))}(undef, Nsim_2)
    Xp = @views X[1:Nsim_2, :]
    Xm = @views X[(Nsim_2+1):end, :]
    @inbounds for j = 1:Nstep
        tmp_ = incremental_integral(μ, (j - 1) * dt, dt)
        randn!(mcBaseData.parallelMode.rng, Z)
        @views @. Xp[:, j+1] = Xp[:, j] + tmp_ + stddev_bm * Z
        @views @. Xm[:, j+1] = Xm[:, j] + tmp_ - stddev_bm * Z
    end

    nothing
end
