
function add_jump_matrix!(X, mcProcess::finiteActivityProcess, mcBaseData::MonteCarloConfiguration{type1, type2, type3, type4, type5}, T::numb) where {finiteActivityProcess <: FiniteActivityProcess, numb <: Number, type1 <: Number, type2 <: Number, type3 <: AbstractMonteCarloMethod, type4 <: BaseMode, type5 <: Random.AbstractRNG}
    Nsim = mcBaseData.Nsim
    Nstep = mcBaseData.Nstep
    λ = mcProcess.λ
    dt = T / Nstep
    @inbounds for i = 1:Nsim
        t_i = randexp(mcBaseData.rng) / λ
        while t_i < T
            jump_size = compute_jump_size(mcProcess, mcBaseData)
            jump_idx = ceil(UInt32, t_i / dt) + 1
            @views X[i, jump_idx:end] .+= jump_size #add jump component
            t_i += randexp(mcBaseData.rng) / λ
        end
    end
end

function simulate!(X, mcProcess::finiteActivityProcess, rfCurve::AbstractZeroRateCurve, mcBaseData::MonteCarloConfiguration{type1, type2, type3, type4, type5}, T::numb) where {finiteActivityProcess <: FiniteActivityProcess, numb <: Number, type1 <: Number, type2 <: Number, type3 <: AbstractMonteCarloMethod, type4 <: BaseMode, type5 <: Random.AbstractRNG}
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
