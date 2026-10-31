
function simulate!(X, mcProcess::finiteActivityProcess, rfCurve::AbstractZeroRateCurve, mcBaseData::MonteCarloConfiguration{type1, type2, type3, <:AbstractCudaMode}, T::numb) where {finiteActivityProcess <: FiniteActivityProcess, numb <: Number, type1 <: Number, type2 <: Number, type3 <: AbstractMonteCarloMethod}
    r = rfCurve.r
    d = dividend(mcProcess)
    Nsim = mcBaseData.Nsim
    Nstep = mcBaseData.Nstep
    σ = mcProcess.σ
    λ = mcProcess.λ

    ChainRulesCore.@ignore_derivatives @assert T > 0.0
    ## Simulate
    # r-d-psi(-i)
    drift_rn = (r - d) - σ^2 / 2 - λ * compute_drift(mcProcess)
    simulate!(X, BrownianMotion(σ, drift_rn), mcBaseData, T)
    # The jump times are drawn on the host, and writing them to X would launch a kernel for each jump.
    X_incr = zeros(eltype(X), Nsim, Nstep + 1)
    add_jump_matrix!(X_incr, mcProcess, mcBaseData, T)
    ## Conclude
    S0 = mcProcess.underlying.S0
    X_incr_cu = CuArray(X_incr)
    @. X = S0 * exp(X + X_incr_cu)
    return
end
