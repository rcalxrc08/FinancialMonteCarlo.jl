using Sobol, StatsFuns

# What the simulate! methods share, so that one body covers antithetic variates, constant and curve rates, and the
# devices: step_drifts handles the rates, path_sides the antithetic variates, and the arrays come from the mode
# (array_type): a device adds its method of array_type and, for the arrays of its own type, of draw_normals! and
# subordinator_host.

# The paths of X with the sign of their normal draws: all of them or, with antithetic variates, the two halves
path_sides(X, ::MonteCarloConfiguration) = ((X, Int8(1)),)

function path_sides(X, mcBaseData::MonteCarloConfiguration{<:AntitheticMC})
    Nsim_2 = div(mcBaseData.Nsim, 2)
    return ((view(X, 1:Nsim_2, :), Int8(1)), (view(X, (Nsim_2+1):mcBaseData.Nsim, :), Int8(-1)))
end

# The drift of each time step: a constant rate is one number, only a curve needs a vector
step_drifts(μ::Number, dt::Number, ::Integer) = μ * dt
step_drifts(μ::CurveType, dt::Number, Nstep::Integer) = [incremental_integral(μ, (j - 1) * dt, dt) for j = 1:Nstep]
drift_of_step(drifts::Number, ::Integer) = drifts
drift_of_step(drifts::AbstractVector, j::Integer) = drifts[j]

draws_per_side(sides) = size(sides[1][1], 1)
# The array type of a mode: the CUDA extension adds CuArray for its modes
array_type(::AbstractMode) = Array
draw_buffer(mcBaseData::MonteCarloConfiguration, numtype, N::Integer) = array_type(mcBaseData.mode){numtype}(undef, N)
# The CUDA extension draws the normals of a CuArray from CURAND
draw_normals!(Z, mcBaseData::MonteCarloConfiguration) = randn!(host_rng(mcBaseData), Z)
# The laws of Distributions sample on the host: the array to fill there, which on the CPU is dt_s itself
subordinator_host(dt_s) = dt_s

function draw_subordinator!(dt_s, dt_s_host, subordinator_, mcBaseData::MonteCarloConfiguration)
    rand!(host_rng(mcBaseData), subordinator_, dt_s_host)
    dt_s === dt_s_host || copyto!(dt_s, dt_s_host)
    return dt_s
end

# The type of a Sobol sequence holds its dimension, known only at run time, so the paths are drawn behind a
# function barrier: path!(i, normals) fills the path i from the normals of its point of the sequence.
on_sobol_paths(path!, dimension::Integer, mcBaseData::MonteCarloConfiguration, numtype::Type) = sobol_paths!(path!, SobolSeq(dimension), Vector{numtype}(undef, dimension), mcBaseData)

# Sobol.skip(seq, n) skips 2^m - 1 points, 2^m <= n + 1 < 2^(m + 1), by drawing every one of them: with n = Nstep * Nsim
# that was 32 times the points of the simulation, most of its time. After k points the state of the sequence has a closed
# form (the Gray code of k selects the direction numbers, every coordinate has its fixed point after bit floor(log2(k))),
# so the same state is set here in O(dimension * log2(k)), with the same points afterwards.
function sobol_skip!(seq::SobolSeq, n::Integer)
    n <= 0 && return seq
    k = (one(UInt64) << floor(Int, log2(n + 1))) - one(UInt64)
    k == 0 && return seq
    B = 63 - leading_zeros(k)
    gray = k ⊻ (k >> 1)
    for i in eachindex(seq.x)
        x = zero(UInt32)
        for c = 0:B
            if isodd(gray >> c)
                x ⊻= seq.m[i, c+1] << (B - c)
            end
        end
        seq.x[i] = x
        seq.b[i] = B
    end
    seq.n = k
    return seq
end

function sobol_paths!(path!, seq, normals, mcBaseData::MonteCarloConfiguration)
    sobol_skip!(seq, mcBaseData.Nstep * mcBaseData.Nsim)
    for i = 1:mcBaseData.Nsim
        next!(seq, normals)
        @. normals = norminvcdf(normals)
        path!(i, normals)
    end
    return
end
