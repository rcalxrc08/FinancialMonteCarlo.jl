"""
Struct for Brownian Motion

		bmProcess=BrownianMotion(σ::num1,μ::num2) where {num1,num2 <: Number}
	
Where:

		σ	=	volatility of the process.
		μ	=	drift of the process.
"""
struct BrownianMotion{num <: Number, num1 <: Number, numtype <: Number} <: AbstractMonteCarloEngine{numtype}
    σ::num
    μ::num1
    function BrownianMotion(σ::num, μ::num1) where {num <: Number, num1 <: Number}
        @assert σ > 0 "Volatility must be positive"
        zero_typed = zero(num) + zero(num1)
        return new{num, num1, typeof(zero_typed)}(σ, μ)
    end
end

export BrownianMotion;

# function simulate_path!(X, mcProcess::BrownianMotion, mcBaseData::SerialMonteCarloConfig, T::numb) where {numb <: Number}
#     Nstep = mcBaseData.Nstep
#     σ = mcProcess.σ
#     μ = mcProcess.μ
#     @assert T > 0
#     dt = T / Nstep
#     mean_bm = μ * dt
#     stddev_bm = σ * sqrt(dt)
#     isDualZero = mean_bm * stddev_bm * 0
#     @views X[1] = isDualZero
#     @inbounds for j = 1:Nstep
#         @views X[j+1] = X[j] + mean_bm + stddev_bm * randn(mcBaseData.parallelMode.rng)
#     end
#     return
# end

function simulate!(X, mcProcess::BrownianMotion, mcBaseData::SerialMonteCarloConfig, T::numb) where {numb <: Number}
    Nsim = mcBaseData.Nsim
    Nstep = mcBaseData.Nstep
    σ = mcProcess.σ
    μ = mcProcess.μ
    @assert T > 0
    dt = T / Nstep
    mean_bm = μ * dt
    stddev_bm = σ * sqrt(dt)
    isDualZero = mean_bm * stddev_bm * 0
    view(X, :, 1) .= isDualZero
    Z = Array{typeof(get_rng_type(isDualZero))}(undef, Nsim)
    @inbounds for j = 1:Nstep
        randn!(mcBaseData.parallelMode.rng, Z)
        @views @. X[:, j+1] = X[:, j] + mean_bm + stddev_bm * Z
    end
    return
end

function simulate!(X, mcProcess::BrownianMotion, mcBaseData::SerialAntitheticMonteCarloConfig, T::numb) where {numb <: Number}
    Nsim = mcBaseData.Nsim
    Nstep = mcBaseData.Nstep
    σ = mcProcess.σ
    μ = mcProcess.μ
    @assert T > 0
    dt = T / Nstep
    mean_bm = μ * dt
    stddev_bm = σ * sqrt(dt)
    isDualZero = mean_bm * stddev_bm * 0
    view(X, :, 1) .= isDualZero
    Nsim_2 = div(Nsim, 2)
    zero_typed = predict_output_type_zero(σ, μ)
    Z = Array{typeof(get_rng_type(zero_typed))}(undef, Nsim_2)
    Xp = @views X[1:Nsim_2, :]
    Xm = @views X[(Nsim_2+1):end, :]
    @inbounds for j = 1:Nstep
        randn!(mcBaseData.parallelMode.rng, Z)
        @. @views Xp[:, j+1] = Xp[:, j] + mean_bm + stddev_bm * Z
        @. @views Xm[:, j+1] = Xm[:, j] + mean_bm - stddev_bm * Z
    end

    nothing
end
