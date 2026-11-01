using Sobol, StatsFuns

# What the simulate! methods share, so that one body covers antithetic variates, constant and curve rates, and the
# devices: step_drifts handles the rates, path_sides the antithetic variates, and a device only adds its methods of
# draw_buffer, draw_normals! and subordinator_host.
const PseudoRandomConfig = Union{SerialMonteCarloConfig, CudaMonteCarloConfig, CudaAntitheticMonteCarloConfig}

# The paths of X with the sign of their normal draws: all of them or, with antithetic variates, the two halves
path_sides(X, ::MonteCarloConfiguration) = ((X, Int8(1)),)

function path_sides(X, mcBaseData::MonteCarloConfiguration{<:Integer, <:Integer, <:AntitheticMC})
    Nsim_2 = div(mcBaseData.Nsim, 2)
    return ((view(X, 1:Nsim_2, :), Int8(1)), (view(X, (Nsim_2+1):mcBaseData.Nsim, :), Int8(-1)))
end

# The drift of each time step: a constant rate is one number, only a curve needs a vector
step_drifts(μ::Number, dt::Number, ::Integer) = μ * dt
step_drifts(μ::CurveType, dt::Number, Nstep::Integer) = [incremental_integral(μ, (j - 1) * dt, dt) for j = 1:Nstep]
drift_of_step(drifts::Number, ::Integer) = drifts
drift_of_step(drifts::AbstractVector, j::Integer) = drifts[j]

draws_per_side(sides) = size(sides[1][1], 1)
draw_buffer(::MonteCarloConfiguration, numtype, N::Integer) = Array{numtype}(undef, N)
draw_normals!(Z, mcBaseData::MonteCarloConfiguration) = randn!(mcBaseData.parallelMode.rng, Z)
# The laws of Distributions sample on the host: the array to fill there, which on the CPU is dt_s itself
subordinator_host(dt_s) = dt_s

function draw_subordinator!(dt_s, dt_s_host, subordinator_, mcBaseData::MonteCarloConfiguration)
    rand!(mcBaseData.parallelMode.rng, subordinator_, dt_s_host)
    dt_s === dt_s_host || copyto!(dt_s, dt_s_host)
    return dt_s
end

function sobol_sequence(dimension::Integer, mcBaseData::MonteCarloConfiguration)
    seq = SobolSeq(dimension)
    skip(seq, mcBaseData.Nstep * mcBaseData.Nsim)
    return seq
end

function sobol_normals!(vec, seq)
    next!(seq, vec)
    @. vec = norminvcdf(vec)
    return vec
end
