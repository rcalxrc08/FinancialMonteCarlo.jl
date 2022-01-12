
function add_jump_matrix!(X, mcProcess::finiteActivityProcess, mcBaseData::AbstractMonteCarloConfiguration, T::numb) where {finiteActivityProcess <: FiniteActivityProcess, numb <: Number}
    Nsim = mcBaseData.Nsim
    Nstep = mcBaseData.Nstep
    λ = mcProcess.λ
    dt = T / Nstep
    @inbounds for i = 1:Nsim
        t_i = randexp(mcBaseData.parallelMode.rng) / λ
        while t_i < T
            jump_size = compute_jump_size(mcProcess, mcBaseData)
            jump_idx = ceil(UInt32, t_i / dt) + 1
            @views @. X[i, jump_idx:end] += jump_size #add jump component
            t_i += randexp(mcBaseData.parallelMode.rng) / λ
        end
    end
end

function simulate!(X, mcProcess::finiteActivityProcess, rfCurve::AbstractZeroRateCurve, mcBaseData::AbstractMonteCarloConfiguration, T::numb) where {finiteActivityProcess <: FiniteActivityProcess, numb <: Number}
    r = rfCurve.r
    d = dividend(mcProcess)
    σ = mcProcess.σ
    λ = mcProcess.λ

    @assert T > 0.0

    ####Simulation
    ## Simulate
    # r-d-psi(-i)
    drift_RN = (r - d) - σ^2 / 2 - λ * compute_drift(mcProcess)
    simulate!(X, BrownianMotion(σ, drift_RN), mcBaseData, T)
    add_jump_matrix!(X, mcProcess, mcBaseData, T)
    ## Conclude
    S0 = mcProcess.underlying.S0
    @. X = S0 * exp(X)
    nothing
end
