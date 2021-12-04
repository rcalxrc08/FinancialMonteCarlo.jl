
function simulate!(X, mcProcess::HestonProcess, rfCurve::ZeroRateCurve, mcBaseData::MonteCarloConfiguration{type1, type2, type3, SerialMode, type5}, T::numb) where {numb <: Number, type1 <: Integer, type2 <: Integer, type3 <: StandardMC, type5 <: Random.AbstractRNG}
    r = rfCurve.r
    S0 = mcProcess.underlying.S0
    d = dividend(mcProcess)
    Nsim = mcBaseData.Nsim
    Nstep = mcBaseData.Nstep
    σ = mcProcess.σ
    σ₀ = mcProcess.σ₀
    λ1 = mcProcess.λ
    κ = mcProcess.κ
    ρ = mcProcess.ρ
    θ = mcProcess.θ
    @assert T > 0.0

    ## Simulate
    κ_s = κ + λ1
    θ_s = κ * θ / (κ + λ1)
    r_d = r - d
    dt = T / Nstep
    zero_drift = incremental_integral(r_d, dt * 0.0, dt)
    isDualZero = S0 * T * σ₀ * κ * θ * λ1 * σ * ρ * 0.0 * zero_drift
    X[:, 1] .= isDualZero
    v_m = [σ₀^2 + isDualZero for _ = 1:Nsim]
    for j = 1:Nstep
        e1 = randn(mcBaseData.rng, Nsim)
        e2 = randn(mcBaseData.rng, Nsim)
        e2 = e1 * ρ + e2 * sqrt(1 - ρ * ρ)
        tmp_ = incremental_integral(r_d, (j - 1) * dt, dt)
        X[:, j+1] = X[:, j] + (tmp_ .- 0.5 * max.(v_m, 0) * dt) + sqrt.(max.(v_m, 0)) * sqrt(dt) .* e1
        v_m += κ_s .* (θ_s .- max.(v_m, 0)) .* dt + σ .* sqrt.(max.(v_m, 0)) .* sqrt(dt) .* e2
    end
    ## Conclude
    X .= S0 .* exp.(X)
    return
end

function simulate!(X, mcProcess::HestonProcess, rfCurve::ZeroRateCurve, mcBaseData::MonteCarloConfiguration{type1, type2, type3, SerialMode, type5}, T::numb) where {numb <: Number, type1 <: Integer, type2 <: Integer, type3 <: AntitheticMC, type5 <: Random.AbstractRNG}
    r = rfCurve.r
    S0 = mcProcess.underlying.S0
    d = dividend(mcProcess)
    Nsim = mcBaseData.Nsim
    Nstep = mcBaseData.Nstep
    σ = mcProcess.σ
    σ₀ = mcProcess.σ₀
    λ1 = mcProcess.λ
    κ = mcProcess.κ
    ρ = mcProcess.ρ
    θ = mcProcess.θ
    @assert T > 0.0

    ####Simulation
    ## Simulate
    κ_s = κ + λ1
    θ_s = κ * θ / (κ + λ1)

    dt = T / Nstep
    isDualZero = S0 * T * σ₀ * κ * θ * λ1 * σ * ρ * 0.0

    X[:, 1] .= isDualZero
    Nsim_2 = div(Nsim, 2)
    r_d = r - d
    Xp = @views X[1:Nsim_2, :]
    Xm = @views X[(Nsim_2+1):end, :]
    v_mp = [σ₀^2 + isDualZero for _ = 1:Nsim_2]
    v_mm = [σ₀^2 + isDualZero for _ = 1:Nsim_2]
    e1 = Array{typeof(get_rng_type(isDualZero))}(undef, Nsim_2)
    e2 = Array{typeof(get_rng_type(isDualZero))}(undef, Nsim_2)
    for j = 1:Nstep
        randn!(mcBaseData.rng, e1)
        randn!(mcBaseData.rng, e2)
        e2 = e1 .* ρ .+ e2 .* sqrt(1 - ρ * ρ)
        tmp_ = incremental_integral(r_d, (j - 1) * dt, dt)
        @views Xp[:, j+1] = Xp[:, j] + (tmp_ .- 0.5 .* max.(v_mp, 0) .* dt) + sqrt.(max.(v_mp, 0)) .* sqrt(dt) .* e1
        @views Xm[:, j+1] = Xm[:, j] + (tmp_ .- 0.5 .* max.(v_mm, 0) .* dt) - sqrt.(max.(v_mm, 0)) .* sqrt(dt) .* e1
        v_mp += κ_s .* (θ_s .- max.(v_mp, 0)) .* dt + (σ * sqrt(dt)) .* sqrt.(max.(v_mp, 0)) .* e2
        v_mm += κ_s .* (θ_s .- max.(v_mm, 0)) .* dt - (σ * sqrt(dt)) .* sqrt.(max.(v_mm, 0)) .* e2
    end

    ## Conclude
    X .= S0 .* exp.(X)
    return
end
