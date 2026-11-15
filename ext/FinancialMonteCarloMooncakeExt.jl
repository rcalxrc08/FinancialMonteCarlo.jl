module FinancialMonteCarloMooncakeExt

using FinancialMonteCarlo, Mooncake
using Mooncake: @zero_derivative, @is_primitive, MinimalCtx, CoDual, Dual, NoRData, primal, tangent, arrayify, zero_fcodual, zero_dual
using FinancialMonteCarlo: bm_step!, exp_paths!, affine_carry!, add_carry!, heston_step!, heston_log_step, heston_variance_step, lnm_step!, lnm_euler, path_range
using FinancialMonteCarlo: add_jump_sizes!, add_jump_sizes_adjoint!, record_jumps, bm_step_adjoint!, exp_paths_adjoint!, affine_carry_adjoint!, add_carry_adjoint!, heston_step_adjoint!, lnm_step_adjoint!, lnm_variance_and_partials
using MuladdMacro: @muladd

# Reseeding the rng ends in a foreigncall for the MersenneTwister (the default rng), and the seed carries no derivative.
@zero_derivative MinimalCtx Tuple{typeof(FinancialMonteCarlo.set_seed!), Any}
# The minimum and maximum of every path reach the barrier payoffs only through comparisons
@zero_derivative MinimalCtx Tuple{typeof(path_range), Any}
# The jumps as random numbers and grid points: their sizes are computed after, from the parameters
@zero_derivative MinimalCtx Tuple{typeof(record_jumps), Any, Any, Any}

# Rules for the kernels of the paths, on arrays of floats (paths.jl, heston.jl, log_normal_mixture.jl, fa_levy.jl).
# Without them Mooncake differentiates their broadcasts element by element, pushing the data of every element on its
# stack; with them each step is one pass over the paths, with the derivatives in closed form (src/models/adjoints.jl).
#
# Anatomy of a rule, for a kernel f!(Y, X, p) that writes Y from the arrays X and the scalars p:
# - `@is_primitive` tells Mooncake to call the rules below instead of differentiating f! itself. The signature restricts
#   them to arrays of IEEE floats: dual numbers, GPU arrays and other types go through the generic code.
# - Forward mode, `frule!!`: every argument is a `Dual` (its value and its tangent, the directional derivative carried
#   along). The rule writes the tangent of the output, dY = ∂f/∂X dX + ∂f/∂p dp, into the tangent of Y, then the value of
#   Y. Tangent first: the formulas read the old values of the arrays the kernel overwrites.
# - Reverse mode, `rrule!!`: every argument is a `CoDual` (its value and its fdata, for an array the buffer in which its
#   adjoint is accumulated). The rule runs the kernel and returns a pullback, which Mooncake calls in the reverse pass
#   with the adjoint of the output already in the fdata of Y. The pullback calls the adjoint of adjoints.jl (it adds the
#   adjoints of the arrays into their fdata, in place, and returns those of the scalars), restores what the kernel
#   overwrote, and returns the rdata: `NoRData()` for the function and for every array (their adjoints are already in
#   their fdata), the adjoint of each scalar argument.
# - Mooncake requires the primal state to be restored in the pullback, as its own in-place rules do: the rule copies,
#   before running the kernel, what the kernel overwrites (Y) and what the adjoint needs that may change afterwards (the
#   normals, which the next step redraws into the same buffer).
const IEEEFloat = Union{Float16, Float32, Float64}
const FloatVector{P} = Union{Vector{P}, SubArray{P, 1, <:Array{P}}}

# bm_step!: dest = src + drift + s Z (one step of a Brownian motion). Forward: dy = dx + ddrift + ds Z + s dZ.
# Reverse: the new column (overwritten, restored) and the normals (redrawn at the next step) are saved.
@is_primitive MinimalCtx Tuple{typeof(bm_step!), FloatVector{P}, FloatVector{P}, P, P, Vector{P}} where {P <: IEEEFloat}
@muladd function Mooncake.frule!!(::Dual{typeof(bm_step!)}, dest::Dual{<:FloatVector{P}}, src::Dual{<:FloatVector{P}}, drift::Dual{P}, s::Dual{P}, Z::Dual{Vector{P}}) where {P <: IEEEFloat}
    y, dy = arrayify(dest)
    x, dx = arrayify(src)
    z, dz = arrayify(Z)
    σ, dσ = primal(s), tangent(s)
    ddrift = tangent(drift)
    @. dy = dx + ddrift + dσ * z + σ * dz
    bm_step!(y, x, primal(drift), σ, z)
    return zero_dual(nothing)
end
@muladd function Mooncake.rrule!!(::CoDual{typeof(bm_step!)}, dest::CoDual{<:FloatVector{P}}, src::CoDual{<:FloatVector{P}}, drift::CoDual{P}, s::CoDual{P}, Z::CoDual{Vector{P}}) where {P <: IEEEFloat}
    y, dy = arrayify(dest)
    x, dx = arrayify(src)
    z, dz = arrayify(Z)
    σ = primal(s)
    y_old, z_saved = copy(y), copy(z)
    bm_step!(y, x, primal(drift), σ, z)
    function bm_step_pb!!(::NoRData)
        ddrift, dσ = bm_step_adjoint!(dx, dz, dy, z_saved, σ)
        y .= y_old
        return NoRData(), NoRData(), NoRData(), ddrift, dσ, NoRData()
    end
    return zero_fcodual(nothing), bm_step_pb!!
end

# exp_paths!: X = S0 exp(X), in place. Forward: dY = (dS0 + S0 dX) exp(X). Reverse: the log paths are saved; the
# adjoint needs the new values (still in X when the pullback runs) and the old ones only if S0 = 0.
@is_primitive MinimalCtx Tuple{typeof(exp_paths!), Array{P}, P} where {P <: IEEEFloat}
@muladd function Mooncake.frule!!(::Dual{typeof(exp_paths!)}, X::Dual{<:Array{P}}, S0::Dual{P}) where {P <: IEEEFloat}
    x, dx = arrayify(X)
    s0, ds0 = primal(S0), tangent(S0)
    @. x = exp(x)
    @. dx = (ds0 + s0 * dx) * x
    x .*= s0
    return zero_dual(nothing)
end
function Mooncake.rrule!!(::CoDual{typeof(exp_paths!)}, X::CoDual{<:Array{P}}, S0::CoDual{P}) where {P <: IEEEFloat}
    x, dx = arrayify(X)
    s0 = primal(S0)
    x_old = copy(x)
    exp_paths!(x, s0)
    function exp_paths_pb!!(::NoRData)
        ds0 = exp_paths_adjoint!(dx, x, x_old, s0)
        x .= x_old
        return NoRData(), NoRData(), ds0
    end
    return zero_fcodual(nothing), exp_paths_pb!!
end

# affine_carry!: X = (a + X) c by column (Bachelier), in place. Forward: dY_ij = (a + X_ij) dc_j + (da + dX_ij) c_j.
# Reverse: the old X is saved (the adjoint of c reads it, and it is restored).
@is_primitive MinimalCtx Tuple{typeof(affine_carry!), Matrix{P}, P, Vector{P}} where {P <: IEEEFloat}
@muladd function Mooncake.frule!!(::Dual{typeof(affine_carry!)}, X::Dual{Matrix{P}}, a::Dual{P}, c::Dual{Vector{P}}) where {P <: IEEEFloat}
    x, dx = arrayify(X)
    cc, dc = arrayify(c)
    a_, da = primal(a), tangent(a)
    @inbounds for j in axes(x, 2), i in axes(x, 1)
        dx[i, j] = (a_ + x[i, j]) * dc[j] + (da + dx[i, j]) * cc[j]
        x[i, j] = (a_ + x[i, j]) * cc[j]
    end
    return zero_dual(nothing)
end
@muladd function Mooncake.rrule!!(::CoDual{typeof(affine_carry!)}, X::CoDual{Matrix{P}}, a::CoDual{P}, c::CoDual{Vector{P}}) where {P <: IEEEFloat}
    x, dx = arrayify(X)
    cc, dc = arrayify(c)
    a_ = primal(a)
    x_old = copy(x)
    affine_carry!(x, a_, cc)
    function affine_carry_pb!!(::NoRData)
        ȧ = affine_carry_adjoint!(dx, dc, x_old, a_, cc)
        x .= x_old
        return NoRData(), NoRData(), ȧ, NoRData()
    end
    return zero_fcodual(nothing), affine_carry_pb!!
end

# add_carry!: X = X + b c by column (the shifted models), in place. Forward: dY_ij = dX_ij + b dc_j + db c_j.
# Reverse: the old X is saved only to be restored (the adjoint does not read it).
@is_primitive MinimalCtx Tuple{typeof(add_carry!), Matrix{P}, P, Vector{P}} where {P <: IEEEFloat}
@muladd function Mooncake.frule!!(::Dual{typeof(add_carry!)}, X::Dual{Matrix{P}}, b::Dual{P}, c::Dual{Vector{P}}) where {P <: IEEEFloat}
    x, dx = arrayify(X)
    cc, dc = arrayify(c)
    b_, db = primal(b), tangent(b)
    @inbounds for j in axes(x, 2), i in axes(x, 1)
        dx[i, j] = dx[i, j] + b_ * dc[j] + db * cc[j]
        x[i, j] = x[i, j] + b_ * cc[j]
    end
    return zero_dual(nothing)
end
@muladd function Mooncake.rrule!!(::CoDual{typeof(add_carry!)}, X::CoDual{Matrix{P}}, b::CoDual{P}, c::CoDual{Vector{P}}) where {P <: IEEEFloat}
    x, dx = arrayify(X)
    cc, dc = arrayify(c)
    b_ = primal(b)
    x_old = copy(x)
    add_carry!(x, b_, cc)
    function add_carry_pb!!(::NoRData)
        ḃ = add_carry_adjoint!(dx, dc, b_, cc)
        x .= x_old
        return NoRData(), NoRData(), ḃ, NoRData()
    end
    return zero_fcodual(nothing), add_carry_pb!!
end

# heston_step!: X_to = X_from + drift - v dt / 2 + √v a e1, v = max(v + κ (θ - v) dt + σ √v a e2, v_floor) in place,
# e2 = ρ e1 + ρ_orthogonal e2_rho; the nine scalars are drift, dt, a, κ, θ, σ, ρ, ρ_orthogonal, v_floor.
# Forward: one loop for the values and the tangents, element by element: the tangent of the log price and of the raw
# variance from the tangents of every input (d√v = dv / (2√v)), then the floor picks the tangent of the raw variance or
# of the floor, as max picks the value. Reverse: the new column, the variance before the step (read by the adjoint and
# restored) and the two normals are saved.
@is_primitive MinimalCtx Tuple{typeof(heston_step!), FloatVector{P}, FloatVector{P}, Vector{P}, Vector{P}, Vector{P}, P, P, P, P, P, P, P, P, P} where {P <: IEEEFloat}
@muladd function Mooncake.frule!!(::Dual{typeof(heston_step!)}, X_to::Dual{<:FloatVector{P}}, X_from::Dual{<:FloatVector{P}}, v::Dual{Vector{P}}, e1::Dual{Vector{P}}, e2_rho::Dual{Vector{P}}, scalars::Vararg{Dual{P}, 9}) where {P <: IEEEFloat}
    y, dy = arrayify(X_to)
    x, dx = arrayify(X_from)
    w, dw = arrayify(v)
    z1, dz1 = arrayify(e1)
    z2, dz2 = arrayify(e2_rho)
    drift, dt, a, κ, θ, σ, ρ, ρo, fl = map(primal, scalars)
    ddrift, ddt, da, dκ, dθ, dσ, dρ, dρo, dfl = map(tangent, scalars)
    @inbounds for i in eachindex(y, x, w, z1, z2)
        vi, dvi = w[i], dw[i]                     # the variance before the step and its tangent
        sv = sqrt(vi)
        dsv = dvi / (2 * sv)                      # tangent of √v
        e2 = z1[i] * ρ + z2[i] * ρo
        de2 = dz1[i] * ρ + z1[i] * dρ + dz2[i] * ρo + z2[i] * dρo
        dy[i] = dx[i] + ddrift - (dvi * dt + vi * ddt) / 2 + (dsv * a + sv * da) * z1[i] + sv * a * dz1[i]
        y[i] = heston_log_step(x[i], vi, z1[i], drift, dt, a)
        raw = heston_variance_step(vi, z1[i], z2[i], dt, a, κ, θ, σ, ρ, ρo)
        # tangent of raw = v + κ (θ - v) dt + σ √v a e2, term by term
        draw = dvi + (dκ * (θ - vi) + κ * (dθ - dvi)) * dt + κ * (θ - vi) * ddt + (dσ * sv + σ * dsv) * a * e2 + σ * sv * (da * e2 + a * de2)
        alive = fl < raw
        w[i] = alive ? raw : fl
        dw[i] = alive ? draw : dfl
    end
    return zero_dual(nothing)
end
@muladd function Mooncake.rrule!!(::CoDual{typeof(heston_step!)}, X_to::CoDual{<:FloatVector{P}}, X_from::CoDual{<:FloatVector{P}}, v::CoDual{Vector{P}}, e1::CoDual{Vector{P}}, e2_rho::CoDual{Vector{P}}, scalars::Vararg{CoDual{P}, 9}) where {P <: IEEEFloat}
    y, dy = arrayify(X_to)
    x, dx = arrayify(X_from)
    w, dw = arrayify(v)
    z1, dz1 = arrayify(e1)
    z2, dz2 = arrayify(e2_rho)
    drift, dt, a, κ, θ, σ, ρ, ρo, fl = map(primal, scalars)
    y_old, w_old, z1_saved, z2_saved = copy(y), copy(w), copy(z1), copy(z2)
    heston_step!(y, x, w, z1, z2, drift, dt, a, κ, θ, σ, ρ, ρo, fl)
    function heston_step_pb!!(::NoRData)
        scalar_adjoints = heston_step_adjoint!(dx, dw, dz1, dz2, dy, w_old, z1_saved, z2_saved, drift, dt, a, κ, θ, σ, ρ, ρo, fl)
        y .= y_old
        w .= w_old
        return NoRData(), NoRData(), NoRData(), NoRData(), NoRData(), NoRData(), scalar_adjoints...
    end
    return zero_fcodual(nothing), heston_step_pb!!
end

# lnm_step!: dest = x + drift - σ² h / 2 + √(σ² h) z with z = sign Z and the local variance σ²(t, x - carry) (σ2_0 at
# t = 0); the scalars are t, h, drift, carry, σ2_0 and the tuples η and logw (their tangents and adjoints are tuples).
# Forward: the tangent of σ² from its partial derivatives (lnm_variance_and_partials), then the tangent of the step.
# Reverse: the new column (restored) and the normals (redrawn) are saved; the old column x is read in the pullback,
# where it holds the same values (the steps after this one write other columns, and exp_paths! restores it first).
@is_primitive MinimalCtx Tuple{typeof(lnm_step!), FloatVector{P}, FloatVector{P}, Vector{P}, Integer, P, P, P, P, P, NTuple{N, P}, NTuple{N, P}} where {N, P <: IEEEFloat}
@muladd function Mooncake.frule!!(::Dual{typeof(lnm_step!)}, dest::Dual{<:FloatVector{P}}, src::Dual{<:FloatVector{P}}, Z::Dual{Vector{P}}, sign::Dual{<:Integer}, t::Dual{P}, h::Dual{P}, drift::Dual{P}, carry::Dual{P}, σ2_0::Dual{P}, η::Dual{NTuple{N, P}}, logw::Dual{NTuple{N, P}}) where {N, P <: IEEEFloat}
    y, dy = arrayify(dest)
    x, dx = arrayify(src)
    z, dz = arrayify(Z)
    sg = primal(sign)
    t_, h_, μ, c, s0, η_, l_ = primal(t), primal(h), primal(drift), primal(carry), primal(σ2_0), primal(η), primal(logw)
    dt_, dh_, dμ, dc, ds0, dη_, dl_ = tangent(t), tangent(h), tangent(drift), tangent(carry), tangent(σ2_0), tangent(η), tangent(logw)
    @inbounds for i in eachindex(y, x, z)
        σ2, py, pt, pη, pl = lnm_variance_and_partials(x[i], t_, c, s0, η_, l_)
        # tangent of the local variance: of σ2_0 at t = 0, else through ψ = x - carry, t, η and logw
        dσ2 = iszero(t_) ? ds0 : py * (dx[i] - dc) + pt * dt_ + sum(map(*, pη, dη_)) + sum(map(*, pl, dl_))
        zi = sg * z[i]
        sq = sqrt(σ2 * h_)
        q = zi * sq / (2 * σ2) - h_ / 2          # ∂y/∂σ²
        qh = zi * sq / (2 * h_) - σ2 / 2          # ∂y/∂h
        dy[i] = dx[i] + dμ + sg * sq * dz[i] + qh * dh_ + q * dσ2
        y[i] = lnm_euler(x[i], t_, h_, μ, zi, c, s0, η_, l_)
    end
    return zero_dual(nothing)
end
@muladd function Mooncake.rrule!!(::CoDual{typeof(lnm_step!)}, dest::CoDual{<:FloatVector{P}}, src::CoDual{<:FloatVector{P}}, Z::CoDual{Vector{P}}, sign::CoDual{<:Integer}, t::CoDual{P}, h::CoDual{P}, drift::CoDual{P}, carry::CoDual{P}, σ2_0::CoDual{P}, η::CoDual{NTuple{N, P}}, logw::CoDual{NTuple{N, P}}) where {N, P <: IEEEFloat}
    y, dy = arrayify(dest)
    x, dx = arrayify(src)
    z, dz = arrayify(Z)
    sg = primal(sign)
    t_, h_, μ, c, s0, η_, l_ = primal(t), primal(h), primal(drift), primal(carry), primal(σ2_0), primal(η), primal(logw)
    y_old, z_saved = copy(y), copy(z)
    lnm_step!(y, x, z, sg, t_, h_, μ, c, s0, η_, l_)
    function lnm_step_pb!!(::NoRData)
        scalar_adjoints = lnm_step_adjoint!(dx, dz, dy, x, z_saved, sg, t_, h_, μ, c, s0, η_, l_)
        y .= y_old
        return NoRData(), NoRData(), NoRData(), NoRData(), NoRData(), scalar_adjoints...
    end
    return zero_fcodual(nothing), lnm_step_pb!!
end

# add_jump_sizes!: X[path[k], step[k]:end] += sizes[k] (the jumps of Kou and Merton), in place. Linear in X and in the
# sizes: the forward rule applies the same kernel to the tangents. Reverse: nothing is saved, the adjoint reads only the
# indices; the pullback restores X by subtracting the sizes (exact up to rounding, and nothing reads X before the
# pullback of the step that wrote it).
@is_primitive MinimalCtx Tuple{typeof(add_jump_sizes!), Matrix{P}, Vector{Int32}, Vector{Int32}, Vector{P}} where {P <: IEEEFloat}
function Mooncake.frule!!(::Dual{typeof(add_jump_sizes!)}, X::Dual{Matrix{P}}, path::Dual{Vector{Int32}}, step::Dual{Vector{Int32}}, sizes::Dual{Vector{P}}) where {P <: IEEEFloat}
    x, dx = arrayify(X)
    s, ds = arrayify(sizes)
    add_jump_sizes!(x, primal(path), primal(step), s)
    add_jump_sizes!(dx, primal(path), primal(step), ds)
    return zero_dual(nothing)
end
function Mooncake.rrule!!(::CoDual{typeof(add_jump_sizes!)}, X::CoDual{Matrix{P}}, path::CoDual{Vector{Int32}}, step::CoDual{Vector{Int32}}, sizes::CoDual{Vector{P}}) where {P <: IEEEFloat}
    x, dx = arrayify(X)
    s, ds = arrayify(sizes)
    p, st = primal(path), primal(step)
    add_jump_sizes!(x, p, st, s)
    function add_jump_sizes_pb!!(::NoRData)
        add_jump_sizes_adjoint!(ds, dx, p, st)
        add_jump_sizes!(x, p, st, -s)
        return NoRData(), NoRData(), NoRData(), NoRData(), NoRData()
    end
    return zero_fcodual(nothing), add_jump_sizes_pb!!
end

end
