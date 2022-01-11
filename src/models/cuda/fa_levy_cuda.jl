
function simulate!(X, mcProcess::finiteActivityProcess, rfCurve::AbstractZeroRateCurve, mcBaseData::MonteCarloConfiguration{type1, type2, type3, CudaMode, type5}, T::numb) where {finiteActivityProcess <: FiniteActivityProcess, numb <: Number, type1 <: Number, type2 <: Number, type3 <: AbstractMonteCarloMethod, type5 <: Random.AbstractRNG}
    r = rfCurve.r
    d = dividend(mcProcess)
    Nsim = mcBaseData.Nsim
    Nstep = mcBaseData.Nstep
    σ = mcProcess.σ
    λ = mcProcess.λ

    @assert T > 0.0
    ## Simulate
    # r-d-psi(-i)
    drift_RN = (r - d) - σ^2 / 2 - λ * compute_drift(mcProcess)
    simulate!(X, BrownianMotion(σ, drift_RN), mcBaseData, T)
    X_incr = zeros(Nsim, Nstep + 1)
    add_jump_matrix!(X_incr, mcProcess, mcBaseData, T)
    ## Conclude
    S0 = mcProcess.underlying.S0
    X_incr_cu = cu(X_incr)
    @. X = S0 * exp(X + X_incr_cu)
end
