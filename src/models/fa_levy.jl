# The jumps of a finite activity Lévy process, added to the log paths in X. The CUDA extension adds a method for a
# CuMatrix, which collects the jumps on the host in a JumpList and adds them with one kernel: the methods for the
# samplings are those of add_jumps!, so they do not clash.
add_jump_matrix!(X, mcProcess::FiniteActivityProcess, mcBaseData::AbstractMonteCarloConfiguration, T::Number) = add_jumps!(X, mcProcess, mcBaseData, T)

# The jumps go to a sink through add_jump!(sink, i, jump_idx, jump_size): the jump of the path i in the step jump_idx - 1,
# which moves the grid points from jump_idx on. The sink is the matrix of the log paths, or a JumpList.
function add_jump!(X::AbstractMatrix, i::Integer, jump_idx::Integer, jump_size)
    @inbounds for j = jump_idx:size(X, 2)
        X[i, j] += jump_size
    end
    return
end

# The jumps of all the paths as a list (path, first grid point, size), to add them where the paths are, e.g. on a GPU
# without a dense matrix of jumps: jumps_by_path sorts them by path, add_jumps!(X, jumps) adds them.
struct JumpList{T}
    path::Vector{Int32}
    step::Vector{Int32}
    size::Vector{T}
end
# n: the expected number of jumps, to allocate the list once
JumpList{T}(n::Integer = 0) where {T} = JumpList{T}(sizehint!(Int32[], n), sizehint!(Int32[], n), sizehint!(T[], n))
function add_jump!(J::JumpList, i::Integer, jump_idx::Integer, jump_size)
    push!(J.path, i)
    push!(J.step, jump_idx)
    push!(J.size, jump_size)
    return
end
# The jumps sorted by path (stably, by counting): those of the path i are offsets[i]:(offsets[i + 1] - 1).
function jumps_by_path(J::JumpList, Nsim::Integer)
    offsets = zeros(Int32, Nsim + 1)
    offsets[1] = 1
    for i in J.path
        offsets[i+1] += 1
    end
    cumsum!(offsets, offsets)
    next = offsets[1:Nsim]
    order = similar(J.path)
    for (k, i) in enumerate(J.path)
        order[next[i]] = k
        next[i] += 1
    end
    return offsets, J.step[order], J.size[order]
end
function add_jumps!(X::AbstractMatrix, J::JumpList)
    for k in eachindex(J.path)
        add_jump!(X, J.path[k], J.step[k], J.size[k])
    end
    return X
end

# The grid point from which a jump at the time t of [0, 1] (as a fraction of the maturity) moves the path: a jump in step
# k moves the grid points from k + 1 on.
jump_index(t, Nstep::Integer) = max(ceil(Int, t * Nstep), 1) + 1

# Plain Monte Carlo: the jump times as sums of exponential spacings, the sizes from compute_jump_size, both from the rng.
function add_jumps!(X, mcProcess::FiniteActivityProcess, mcBaseData::AbstractMonteCarloConfiguration, T::Number)
    Nsim = mcBaseData.Nsim
    Nstep = mcBaseData.Nstep
    λ = mcProcess.λ
    dt = T / Nstep
    @inbounds for i = 1:Nsim
        t_i = draw_exponential(mcBaseData) / λ
        while t_i < T
            add_jump!(X, i, ceil(Int, t_i / dt) + 1, compute_jump_size(mcProcess, mcBaseData))
            t_i += draw_exponential(mcBaseData) / λ
        end
    end
    return
end

# AntitheticMC and SobolMC draw the jumps of a path from uniforms by inversion: the count N_T is the quantile of its
# Poisson law at one uniform, then each jump takes a time (given N_T, the times are independent uniforms) and a size (the
# quantile of the jump law, jump_quantile). A Sobol point has a fixed dimension, so a path has M jumps of uniforms, M
# the Poisson quantile at 1 - sqrt(eps): a path with more jumps draws their times and sizes from the rng, as plain Monte
# Carlo, and a path takes 1 + 2M coordinates after those of its diffusion.
# The Poisson distribution function is tabulated up to M in the type of λT (possibly a dual number: comparing it with a
# uniform compares the values, the count is an integer); logp and logλT continue it above M.
struct JumpUniforms{T <: Number}
    M::Int
    cdf::Vector{T}
    logp::T
    logλT::T
end

"""
    JumpUniforms(λT, U)

The table to draw the number of jumps of a compound Poisson process with mean `λT` from a uniform of type `U`, by
inversion: [`FinancialMonteCarlo.jump_count`](@ref) gives the count. `J.M` is the Poisson quantile at `1 - sqrt(eps(U))`:
with `SobolMC` a path takes `1 + 2 J.M` coordinates (the count, then a time and a size per jump) and draws its jumps
from the rng when there are more. Used by `AntitheticMC` and `SobolMC` for `KouProcess` and `MertonProcess`.
"""
function JumpUniforms(λT::Number, ::Type{U}) where {U <: AbstractFloat}
    tol = sqrt(eps(U))
    logλT = log(λT)
    logp = -λT
    c = exp(logp)
    cdf = [c]
    k = 0
    while c < 1 - tol
        logp_next = logp + logλT - log((k + 1) * one(λT))
        p = exp(logp_next)
        # the distribution function cannot reach 1 - tol in floating point: stop once the terms are negligible
        k + 1 > λT && c + p == c && break
        k += 1
        logp = logp_next
        c += p
        push!(cdf, c)
    end
    return JumpUniforms(length(cdf) - 1, cdf, logp, logλT)
end
jump_uniforms(J::JumpUniforms) = 1 + 2 * J.M

"""
    jump_count(J::JumpUniforms, u)

The number of jumps for the uniform `u`: the smallest `n` with ``P(N_T ≤ n) ≥ u``. Above `J.M` the distribution
function is computed term by term.
"""
function jump_count(J::JumpUniforms, u)
    cdf = J.cdf
    @inbounds for n in eachindex(cdf)
        u <= cdf[n] && return n - 1
    end
    n = J.M
    logp = J.logp
    c = last(cdf)
    while true
        n += 1
        logp += J.logλT - log(n * one(logp))
        c_new = c + exp(logp)
        (u <= c_new || c_new == c) && return n
        c = c_new
    end
end

"""
    jump_quantile(model, u)

The jump size of a `MertonProcess` (normal) or of a `KouProcess` (double exponential: downwards with probability
`1 - p` and rate `λ₋`, upwards with rate `λ₊`) at the uniform `u`, by inversion. Kou's quantile is continuous at
`u = 1 - p` and in `p`. The methods are defined with the models.
"""
function jump_quantile end
# A uniform of a Sobol point or of a reflection can be 0 or 1, where the quantiles are infinite.
open_unit(u) = clamp(u, eps(typeof(u)), 1 - eps(typeof(u)))

# The jumps of one path of a Sobol point, from its uniforms u: the count at u[1], then a time and a size for each jump,
# u[2k] and u[2k + 1]. Above M jumps, the times and sizes come from the rng.
function add_jumps_from_uniforms!(X, i::Integer, u, J::JumpUniforms, mcProcess::FiniteActivityProcess, mcBaseData::AbstractMonteCarloConfiguration)
    Nstep = mcBaseData.Nstep
    N = jump_count(J, open_unit(u[1]))
    if N > J.M
        for _ = 1:N
            add_jump!(X, i, jump_index(draw_uniform(mcBaseData), Nstep), compute_jump_size(mcProcess, mcBaseData))
        end
        return
    end
    @inbounds for k = 1:N
        add_jump!(X, i, jump_index(u[2k], Nstep), jump_quantile(mcProcess, open_unit(u[2k+1])))
    end
    return
end

# Antithetic pairs as the diffusion: the path i + Nsim / 2 takes the uniforms 1 - u of the jumps of the path i. The count
# of the pair comes from one uniform, then each jump from two, drawn only as many as the path with more jumps needs: no
# fixed dimension here, the table is only used for the counts.
function add_jumps!(X, mcProcess::FiniteActivityProcess, mcBaseData::MonteCarloConfiguration{<:AntitheticMC}, T::Number)
    J = JumpUniforms(mcProcess.λ * T, typeof(draw_uniform(mcBaseData)))
    half = div(mcBaseData.Nsim, 2)
    Nstep = mcBaseData.Nstep
    for i = 1:half
        u = open_unit(draw_uniform(mcBaseData))
        N = jump_count(J, u)
        N_anti = jump_count(J, 1 - u)
        for k = 1:max(N, N_anti)
            t = draw_uniform(mcBaseData)
            v = open_unit(draw_uniform(mcBaseData))
            k <= N && add_jump!(X, i, jump_index(t, Nstep), jump_quantile(mcProcess, v))
            k <= N_anti && add_jump!(X, i + half, jump_index(1 - t, Nstep), jump_quantile(mcProcess, 1 - v))
        end
    end
    return
end

# The risk neutral drift of the log price: r - d - ψ(-i)
fa_drift(mcProcess::FiniteActivityProcess, rfCurve::AbstractZeroRateCurve) = (rfCurve.r - dividend(mcProcess)) - mcProcess.σ^2 / 2 - mcProcess.λ * compute_drift(mcProcess)

function simulate!(X, mcProcess::finiteActivityProcess, rfCurve::AbstractZeroRateCurve, mcBaseData::AbstractMonteCarloConfiguration, T::numb) where {finiteActivityProcess <: FiniteActivityProcess, numb <: Number}
    ChainRulesCore.@ignore_derivatives @assert T > 0
    simulate!(X, BrownianMotion(mcProcess.σ, fa_drift(mcProcess, rfCurve)), mcBaseData, T)
    add_jump_matrix!(X, mcProcess, mcBaseData, T)
    S0 = mcProcess.underlying.S0
    @. X = S0 * exp(X)
    nothing
end

# Quasi Monte Carlo: each path is one point of dimension Nstep + 1 + 2M, the first Nstep coordinates give the normals of
# the diffusion (as BrownianMotion), the others the uniforms of the jumps (add_jumps_from_uniforms!).
function simulate!(X, mcProcess::FiniteActivityProcess, rfCurve::AbstractZeroRateCurve, mcBaseData::SobolConfig, T::Number)
    ChainRulesCore.@ignore_derivatives @assert T > 0
    Nstep = mcBaseData.Nstep
    dt = T / Nstep
    stddev_bm = mcProcess.σ * sqrt(dt)
    drifts = step_drifts(fa_drift(mcProcess, rfCurve), dt, Nstep)
    isDualZero = stddev_bm * 0 * drift_of_step(drifts, 1)
    fill!(view(X, :, 1), isDualZero)
    U = typeof(get_rng_type(isDualZero))
    J = JumpUniforms(mcProcess.λ * T, U)
    dimension = Nstep + jump_uniforms(J)
    on_sobol_paths(dimension, mcBaseData, U; normals = Nstep) do i, point
        @inbounds for j = 1:Nstep
            X[i, j+1] = X[i, j] + drift_of_step(drifts, j) + stddev_bm * point[j]
        end
        add_jumps_from_uniforms!(X, i, view(point, (Nstep+1):dimension), J, mcProcess, mcBaseData)
    end
    S0 = mcProcess.underlying.S0
    @. X = S0 * exp(X)
    return
end
