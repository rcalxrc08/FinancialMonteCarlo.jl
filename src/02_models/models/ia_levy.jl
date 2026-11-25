
function simulate!(X, mcProcess::infiniteActivityProcess, und::AbstractUnderlying, rfCurve::AbstractZeroRateCurve, mcBaseData::AbstractMonteCarloConfiguration, T::numb) where {infiniteActivityProcess <: InfiniteActivityProcess, numb <: Number}
    r = rfCurve.r
    d = dividend(und)
    σ = mcProcess.σ
    θ = mcProcess.θ
    Nstep = mcBaseData.Nstep

    check_maturity(T)

    dt = T / Nstep
    ####Simulation
    ## Simulate
    # r-d-psi(-i)
    drift_rn = r - d - compute_drift(mcProcess)
    simulate!(X, SubordinatedBrownianMotion(σ, drift_rn, θ, compute_subordinator_law(mcProcess, dt)), mcBaseData, T)
    ## Conclude
    S0 = und.S0
    exp_paths!(X, S0)
    nothing
end
