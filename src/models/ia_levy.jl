
function simulate!(X, mcProcess::infiniteActivityProcess, rfCurve::AbstractZeroRateCurve, mcBaseData::AbstractMonteCarloConfiguration, T::numb) where {infiniteActivityProcess <: InfiniteActivityProcess, numb <: Number}
    r = rfCurve.r
    d = dividend(mcProcess)
    σ = mcProcess.σ
    θ = mcProcess.θ
    Nstep = mcBaseData.Nstep

    ChainRulesCore.@ignore_derivatives check_domain(T > 0, T, "Time to Maturity must be positive")

    dt = T / Nstep
    ####Simulation
    ## Simulate
    # r-d-psi(-i)
    drift_rn = r - d - compute_drift(mcProcess)
    simulate!(X, SubordinatedBrownianMotion(σ, drift_rn, θ, compute_subordinator_law(mcProcess, dt)), mcBaseData, T)
    ## Conclude
    S0 = mcProcess.underlying.S0
    exp_paths!(X, S0)
    nothing
end
