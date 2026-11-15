module FinancialMonteCarloEnzymeCoreExt

using FinancialMonteCarlo
using EnzymeCore: EnzymeRules, Annotation, Const, Active, Duplicated
using EnzymeCore.EnzymeRules: RevConfigWidth, AugmentedReturn, overwritten
using FinancialMonteCarlo: bm_step!, exp_paths!, affine_carry!, add_carry!, heston_step!, lnm_step!
using FinancialMonteCarlo: add_jump_sizes!, add_jump_sizes_adjoint!, bm_step_adjoint!, exp_paths_adjoint!, affine_carry_adjoint!, add_carry_adjoint!, heston_step_adjoint!, lnm_step_adjoint!

# The seed and the random draws carry no derivative. Reseeding a generator can end in a foreigncall, and Enzyme cannot
# analyse the scalar draws of the MersenneTwister (its state is read as floats and as integers), so they are inactive:
# they run in the forward pass and their results are constants.
EnzymeRules.inactive(::typeof(FinancialMonteCarlo.set_seed!), args...) = nothing
EnzymeRules.inactive(::typeof(FinancialMonteCarlo.draw_uniform), args...) = nothing
EnzymeRules.inactive(::typeof(FinancialMonteCarlo.draw_exponential), args...) = nothing
EnzymeRules.inactive(::typeof(FinancialMonteCarlo.draw_normal), args...) = nothing
# The minimum and maximum of every path reach the barrier payoffs only through comparisons
EnzymeRules.inactive(::typeof(FinancialMonteCarlo.path_range), args...) = nothing
# The jumps as random numbers and grid points: their sizes are computed after, from the parameters
EnzymeRules.inactive(::typeof(FinancialMonteCarlo.record_jumps), args...) = nothing
# and the record itself never carries derivatives
EnzymeRules.inactive_type(::Type{<:FinancialMonteCarlo.JumpRecord}) = true

# Reverse rules for the kernels of the paths, on arrays of floats: the adjoints in closed form of adjoints.jl, instead of
# Enzyme's caching of every value the broadcasts overwrite. Forward mode has no rule: Enzyme differentiates the kernels
# themselves, which is already fast.
#
# Anatomy of a rule, for a kernel f!(Y, X, p) that writes Y from the arrays X and the scalars p:
# - The arguments come annotated: an array is `Duplicated` (its value `.val` and its shadow `.dval`, the buffer of its
#   adjoint) or `Const` (no derivative); a scalar is `Active` (its adjoint is returned by `reverse`) or `Const`.
# - `augmented_primal` runs the kernel in the forward pass and returns a tape: the values the adjoint will need. Enzyme
#   does not restore anything for a custom rule, so what may change before the reverse pass is copied: what the kernel
#   overwrites itself (the variance of Heston), and the arguments that `overwritten(config)` reports as overwritten
#   later by the caller (the normals, redrawn into the same buffer at the next step). The others are kept by reference.
# - `reverse` runs in the reverse pass, with the adjoint of Y in its shadow. It calls the adjoint of adjoints.jl, which
#   adds the adjoints of the arrays into their shadows in place, and returns one entry per argument: the adjoint of each
#   Active scalar, `nothing` for the others. A Const output (Y) means no derivative flows: nothing to do.
const IEEEFloat = Union{Float16, Float32, Float64}
const FloatVector{P} = Union{Vector{P}, SubArray{P, 1, <:Array{P}}}

shadow(x::Duplicated) = x.dval
shadow(::Const) = nothing
# the derivative with respect to an active scalar, nothing for a constant one
scalar(::Active, v) = v
scalar(::Const, v) = nothing
# A copy of the value of an argument if the caller overwrites it after the call. On Julia 1.10 the overwritten flags of
# Enzyme miss some later writes (the column of the log paths read by lnm_step!, which exp_paths! overwrites), with the
# same Enzyme and EnzymeCore that get them right on 1.11 and later: there everything is copied (to drop with 1.10, #66).
saved(config, k, x) = (VERSION < v"1.11" || overwritten(config)[k]) ? copy(x) : x
const ArrayArg{T} = Union{Const{<:T}, Duplicated{<:T}}
const ScalarArg{P} = Union{Const{P}, Active{P}}

# bm_step!: dest = src + drift + s Z. Tape: the normals (only needed for the adjoint of s), copied if redrawn later.
function EnzymeRules.augmented_primal(config::RevConfigWidth{1}, ::Const{typeof(bm_step!)}, ::Type{<:Const}, dest::ArrayArg{FloatVector{P}}, src::ArrayArg{FloatVector{P}}, drift::ScalarArg{P}, s::ScalarArg{P}, Z::ArrayArg{Vector{P}}) where {P <: IEEEFloat}
    tape = s isa Active ? saved(config, 6, Z.val) : nothing
    bm_step!(dest.val, src.val, drift.val, s.val, Z.val)
    return AugmentedReturn(nothing, nothing, tape)
end
function EnzymeRules.reverse(config::RevConfigWidth{1}, ::Const{typeof(bm_step!)}, ::Type{<:Const}, z, dest::ArrayArg{FloatVector{P}}, src::ArrayArg{FloatVector{P}}, drift::ScalarArg{P}, s::ScalarArg{P}, Z::ArrayArg{Vector{P}}) where {P <: IEEEFloat}
    dest isa Const && return (nothing, nothing, nothing, nothing, nothing)
    ddrift, ds = bm_step_adjoint!(shadow(src), shadow(Z), dest.dval, z === nothing ? Z.val : z, s.val)
    return (nothing, nothing, scalar(drift, ddrift), scalar(s, ds), nothing)
end

# exp_paths!: X = S0 exp(X), in place. Tape: the new values (the adjoint multiplies by them), copied only if X is
# overwritten later, and the old values only if S0 = 0.
function EnzymeRules.augmented_primal(config::RevConfigWidth{1}, ::Const{typeof(exp_paths!)}, ::Type{<:Const}, X::ArrayArg{Array{P}}, S0::ScalarArg{P}) where {P <: IEEEFloat}
    x_old = iszero(S0.val) ? copy(X.val) : nothing
    exp_paths!(X.val, S0.val)
    return AugmentedReturn(nothing, nothing, (x_old, saved(config, 2, X.val)))
end
function EnzymeRules.reverse(config::RevConfigWidth{1}, ::Const{typeof(exp_paths!)}, ::Type{<:Const}, tape, X::ArrayArg{Array{P}}, S0::ScalarArg{P}) where {P <: IEEEFloat}
    X isa Const && return (nothing, nothing)
    x_old, x_new = tape
    ds0 = exp_paths_adjoint!(X.dval, x_new, x_old, S0.val)
    return (nothing, scalar(S0, ds0))
end

# affine_carry!: X = (a + X) c by column (Bachelier), in place. Tape: the old X (the adjoint of c reads it) and c.
function EnzymeRules.augmented_primal(config::RevConfigWidth{1}, ::Const{typeof(affine_carry!)}, ::Type{<:Const}, X::ArrayArg{Matrix{P}}, a::ScalarArg{P}, c::ArrayArg{Vector{P}}) where {P <: IEEEFloat}
    tape = (copy(X.val), saved(config, 4, c.val))
    affine_carry!(X.val, a.val, c.val)
    return AugmentedReturn(nothing, nothing, tape)
end
function EnzymeRules.reverse(config::RevConfigWidth{1}, ::Const{typeof(affine_carry!)}, ::Type{<:Const}, tape, X::ArrayArg{Matrix{P}}, a::ScalarArg{P}, c::ArrayArg{Vector{P}}) where {P <: IEEEFloat}
    X isa Const && return (nothing, nothing, nothing)
    x_old, cc = tape
    da = affine_carry_adjoint!(X.dval, shadow(c), x_old, a.val, cc)
    return (nothing, scalar(a, da), nothing)
end

# add_carry!: X = X + b c by column (the shifted models), in place. Tape: c only.
function EnzymeRules.augmented_primal(config::RevConfigWidth{1}, ::Const{typeof(add_carry!)}, ::Type{<:Const}, X::ArrayArg{Matrix{P}}, b::ScalarArg{P}, c::ArrayArg{Vector{P}}) where {P <: IEEEFloat}
    tape = saved(config, 4, c.val)
    add_carry!(X.val, b.val, c.val)
    return AugmentedReturn(nothing, nothing, tape)
end
function EnzymeRules.reverse(config::RevConfigWidth{1}, ::Const{typeof(add_carry!)}, ::Type{<:Const}, cc, X::ArrayArg{Matrix{P}}, b::ScalarArg{P}, c::ArrayArg{Vector{P}}) where {P <: IEEEFloat}
    X isa Const && return (nothing, nothing, nothing)
    db = add_carry_adjoint!(X.dval, shadow(c), b.val, cc)
    return (nothing, scalar(b, db), nothing)
end

# heston_step!: one Euler step of Heston (X_to written, v updated in place). Tape: the variance before the step (always
# copied, the kernel overwrites it) and the two normals.
function EnzymeRules.augmented_primal(config::RevConfigWidth{1}, ::Const{typeof(heston_step!)}, ::Type{<:Const}, X_to::ArrayArg{FloatVector{P}}, X_from::ArrayArg{FloatVector{P}}, v::ArrayArg{Vector{P}}, e1::ArrayArg{Vector{P}}, e2_rho::ArrayArg{Vector{P}}, scalars::Vararg{ScalarArg{P}, 9}) where {P <: IEEEFloat}
    tape = (copy(v.val), saved(config, 5, e1.val), saved(config, 6, e2_rho.val))
    heston_step!(X_to.val, X_from.val, v.val, e1.val, e2_rho.val, map(x -> x.val, scalars)...)
    return AugmentedReturn(nothing, nothing, tape)
end
function EnzymeRules.reverse(config::RevConfigWidth{1}, ::Const{typeof(heston_step!)}, ::Type{<:Const}, tape, X_to::ArrayArg{FloatVector{P}}, X_from::ArrayArg{FloatVector{P}}, v::ArrayArg{Vector{P}}, e1::ArrayArg{Vector{P}}, e2_rho::ArrayArg{Vector{P}}, scalars::Vararg{ScalarArg{P}, 9}) where {P <: IEEEFloat}
    (X_to isa Const || v isa Const) && return (nothing, nothing, nothing, nothing, nothing, ntuple(_ -> nothing, Val(9))...)
    w_old, z1, z2 = tape
    adjoints = heston_step_adjoint!(shadow(X_from), v.dval, shadow(e1), shadow(e2_rho), X_to.dval, w_old, z1, z2, map(x -> x.val, scalars)...)
    return (nothing, nothing, nothing, nothing, nothing, map(scalar, scalars, adjoints)...)
end

# lnm_step!: one Euler step of the lognormal mixture. Tape: the old column of the log paths (the local variance and its
# partials are recomputed from it in the reverse pass) and the normals; `sign` is an integer, always Const.
function EnzymeRules.augmented_primal(config::RevConfigWidth{1}, ::Const{typeof(lnm_step!)}, ::Type{<:Const}, dest::ArrayArg{FloatVector{P}}, src::ArrayArg{FloatVector{P}}, Z::ArrayArg{Vector{P}}, sign::Const{<:Integer}, t::ScalarArg{P}, h::ScalarArg{P}, drift::ScalarArg{P}, carry::ScalarArg{P}, σ2_0::ScalarArg{P}, η::Union{Const{NTuple{N, P}}, Active{NTuple{N, P}}}, logw::Union{Const{NTuple{N, P}}, Active{NTuple{N, P}}}) where {N, P <: IEEEFloat}
    tape = (saved(config, 3, src.val), saved(config, 4, Z.val))
    lnm_step!(dest.val, src.val, Z.val, sign.val, t.val, h.val, drift.val, carry.val, σ2_0.val, η.val, logw.val)
    return AugmentedReturn(nothing, nothing, tape)
end
function EnzymeRules.reverse(config::RevConfigWidth{1}, ::Const{typeof(lnm_step!)}, ::Type{<:Const}, tape, dest::ArrayArg{FloatVector{P}}, src::ArrayArg{FloatVector{P}}, Z::ArrayArg{Vector{P}}, sign::Const{<:Integer}, t::ScalarArg{P}, h::ScalarArg{P}, drift::ScalarArg{P}, carry::ScalarArg{P}, σ2_0::ScalarArg{P}, η::Union{Const{NTuple{N, P}}, Active{NTuple{N, P}}}, logw::Union{Const{NTuple{N, P}}, Active{NTuple{N, P}}}) where {N, P <: IEEEFloat}
    dest isa Const && return ntuple(_ -> nothing, Val(11))
    x, z = tape
    dt, dh, ddrift, dcarry, dσ2_0, dη, dlogw = lnm_step_adjoint!(shadow(src), shadow(Z), dest.dval, x, z, sign.val, t.val, h.val, drift.val, carry.val, σ2_0.val, η.val, logw.val)
    return (nothing, nothing, nothing, nothing, scalar(t, dt), scalar(h, dh), scalar(drift, ddrift), scalar(carry, dcarry), scalar(σ2_0, dσ2_0), scalar(η, dη), scalar(logw, dlogw))
end

# add_jump_sizes!: X[path[k], step[k]:end] += sizes[k] (the jumps of Kou and Merton). Tape: the indices only (the
# adjoint of a size is the sum of the adjoints of the grid points it moved).
function EnzymeRules.augmented_primal(config::RevConfigWidth{1}, ::Const{typeof(add_jump_sizes!)}, ::Type{<:Const}, X::ArrayArg{Matrix{P}}, path::Const{Vector{Int32}}, step::Const{Vector{Int32}}, sizes::ArrayArg{Vector{P}}) where {P <: IEEEFloat}
    add_jump_sizes!(X.val, path.val, step.val, sizes.val)
    return AugmentedReturn(nothing, nothing, (saved(config, 3, path.val), saved(config, 4, step.val)))
end
function EnzymeRules.reverse(config::RevConfigWidth{1}, ::Const{typeof(add_jump_sizes!)}, ::Type{<:Const}, tape, X::ArrayArg{Matrix{P}}, path::Const{Vector{Int32}}, step::Const{Vector{Int32}}, sizes::ArrayArg{Vector{P}}) where {P <: IEEEFloat}
    (X isa Const || sizes isa Const) && return (nothing, nothing, nothing, nothing)
    add_jump_sizes_adjoint!(sizes.dval, X.dval, tape...)
    return (nothing, nothing, nothing, nothing)
end

end
