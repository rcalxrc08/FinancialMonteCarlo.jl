########### Curve Section

# A rate that depends on time: a Curve, or the difference of two rates (CurveDifference). Models only integrate them
# (integral, incremental_integral, step_drifts), so the same code covers constant and time dependent rates.
abstract type AbstractCurve end

"""
    Curve(r, T)

Curve of instantaneous rates with values `r` at the increasing times `T` (two vectors or two tuples; or `r` on a uniform
grid from 0 to a number `T`), linearly interpolated in between and flat outside. Models only use its integrals over time
intervals, both for the drift and for discounting. Curves combine with `+`, `-`, `*` and `/` with numbers (on the
values); two curves can be added and subtracted, which gives a [`FinancialMonteCarlo.CurveDifference`](@ref) (#81).

The nodes are `SVector`s (StaticArrays), so a curve is an isbits value: a model on a curve can be differentiated by
Enzyme without runtime activity and captured by a GPU kernel (#79). The number of nodes is part of the type: a curve
built from vectors has a type known only at run time, so build it once, outside the code that is differentiated or run
in a loop, or from tuples, whose length is part of their type.
"""
struct Curve{N, num1 <: Number, num2 <: Number} <: AbstractCurve
    times::SVector{N, num1}
    values::SVector{N, num2}
    function Curve{N, num1, num2}(times::SVector{N, num1}, values::SVector{N, num2}) where {N, num1 <: Number, num2 <: Number}
        ChainRulesCore.@ignore_derivatives begin
            N >= 1 || throw(ArgumentError("a curve needs at least one node"))
            issorted(times; lt = <=) || throw(ArgumentError("the times of a curve must be increasing"))
        end
        return new{N, num1, num2}(times, values)
    end
    # the times of an existing curve, already checked (map)
    Curve{N, num1, num2}(times::SVector{N, num1}, values::SVector{N, num2}, ::Val{:unchecked}) where {N, num1 <: Number, num2 <: Number} = new{N, num1, num2}(times, values)
end
curve_from_nodes(times::SVector{N, num1}, values::SVector{N, num2}, checked...) where {N, num1 <: Number, num2 <: Number} = Curve{N, num1, num2}(times, values, checked...)
curve_from_nodes(times::NTuple{N, Number}, values::NTuple{N, Number}) where {N} = curve_from_nodes(SVector(times), SVector(values))
# from vectors the number of nodes is known only at run time
function curve_from_nodes(times::AbstractVector, values::AbstractVector)
    ChainRulesCore.@ignore_derivatives length(times) == length(values) || throw(ArgumentError("a curve needs one value per time"))
    N = length(times)
    return curve_from_nodes(SVector{N}(times), SVector{N}(values))
end

# N times from 0 to T
uniform_times(T::Number, N::Integer) = N == 1 ? [zero(T / 1)] : [(i - 1) * (T / (N - 1)) for i = 1:N]
uniform_times(T::Number, ::Val{N}) where {N} = N == 1 ? SVector(zero(T / 1)) : SVector(ntuple(i -> (i - 1) * (T / (N - 1)), Val(N)))
Curve(r_::AbstractVector{<:Number}, T::Number) = curve_from_nodes(uniform_times(T, length(r_)), r_)
Curve(r_::NTuple{N, Number}, T::Number) where {N} = curve_from_nodes(Tuple(uniform_times(T, Val(N))), r_)
Curve(r_::AbstractVector{<:Number}, T::AbstractVector{<:Number}) = curve_from_nodes(T, r_)
Curve(r_::NTuple{N, Number}, T::NTuple{N, Number}) where {N} = curve_from_nodes(T, r_)

"""
    ImpliedCurve(r, T)

Piecewise linear curve of instantaneous rates, starting from zero at time zero, whose integral up to
each `T[i]` equals `r[i] * T[i]`: use it to turn zero rates (or dividend yields) quoted at maturities
`T` into a [`FinancialMonteCarlo.Curve`](@ref). `r` and `T` are two vectors or two tuples; with a scalar `T` the
maturities are `T/N, 2T/N, ..., T`.
"""
function ImpliedCurve(r_::AbstractVector{num1}, T::AbstractVector{num2}) where {num1 <: Number, num2 <: Number}
    ChainRulesCore.@ignore_derivatives check_dimension(length(r_) == length(T), "One rate per maturity")
    ChainRulesCore.@ignore_derivatives check_argument(length(r_) >= 1, "At least one rate")
    # the integral of the piecewise linear curve up to T_i is r_i T_i: on each interval the mean of the curve, the mean of
    # its two end values, is the increment of the integral over the length of the interval
    times = prepend(zero(num2), T)
    means = diff(prepend(zero(num1) * zero(num2), r_ .* T)) ./ diff(times)
    rates = accumulate((previous, mean_) -> 2 * mean_ - previous, means; init = zero(eltype(means)))
    return curve_from_nodes(times, prepend(zero(eltype(rates)), rates))
end
ImpliedCurve(r_::NTuple{N, Number}, T::NTuple{N, Number}) where {N} = ImpliedCurve(SVector(r_), SVector(T))
ImpliedCurve(r_::AbstractVector{<:Number}, T::Number) = ImpliedCurve(r_, [i * (T / length(r_)) for i = 1:length(r_)])
# x followed by the elements of v: a static vector stays static
prepend(x, v::AbstractVector) = vcat(x, v)
prepend(x, v::SVector) = pushfirst(v, x)

Base.keys(curve::Curve) = curve.times
Base.values(curve::Curve) = curve.values
Base.length(::Curve{N}) where {N} = N
Base.map(f, curve::Curve) = curve_from_nodes(curve.times, map(f, curve.values), Val(:unchecked))
Base.:(==)(x::Curve, y::Curve) = x.times == y.times && x.values == y.values
Base.isequal(x::Curve, y::Curve) = isequal(x.times, y.times) && isequal(x.values, y.values)
Base.hash(curve::Curve, h::UInt) = hash(curve.values, hash(curve.times, hash(:Curve, h)))
Base.show(io::IO, curve::Curve) = print(io, "Curve(", join(("$t => $v" for (t, v) in zip(curve.times, curve.values)), ", "), ")")

# The value of the curve at t: linear between the nodes, flat before the first and after the last
function value_at(curve::Curve, t::Number)
    T = curve.times
    r = curve.values
    t <= T[1] && return r[1] + zero(t)
    t >= T[end] && return r[end] + zero(t)
    i = searchsortedlast(T, t)
    @muladd return r[i] + (t - T[i]) / (T[i+1] - T[i]) * (r[i+1] - r[i])
end
value_at(x::Number, ::Number) = x

# Arithmetic with a number, on the values
for op in (:+, :-, :*, :/)
    @eval Base.$op(x::Curve, y::Number) = map(value -> $op(value, y), x)
    @eval Base.$op(x::Number, y::Curve) = map(value -> $op(x, value), y)
end
Base.:-(curve::Curve) = map(-, curve)

"""
    CurveDifference(a, b)

The rate `a - b`, `a` and `b` numbers or curves: what `-` and `+` between two curves return, e.g. the drift `r - d` of a
model whose rate and dividend are both curves. Lazy: its integral is the difference of the integrals of `a` and `b`, exact
without building a curve on the union of their times, and its length is known from the types (#81). Numbers can be added
to it and subtracted from it; curves are not multiplied or divided by each other (the product of two piecewise linear
curves is not piecewise linear).
"""
struct CurveDifference{A <: Union{Number, AbstractCurve}, B <: Union{Number, AbstractCurve}} <: AbstractCurve
    a::A
    b::B
end
Base.:-(x::AbstractCurve, y::AbstractCurve) = CurveDifference(x, y)
Base.:+(x::AbstractCurve, y::AbstractCurve) = CurveDifference(x, -y)
Base.:-(x::CurveDifference) = CurveDifference(x.b, x.a)
Base.:-(x::CurveDifference, y::Number) = CurveDifference(x, y)
Base.:-(x::Number, y::CurveDifference) = CurveDifference(x, y)
Base.:+(x::CurveDifference, y::Number) = CurveDifference(x, -y)
Base.:+(x::Number, y::CurveDifference) = CurveDifference(y, -x)
value_at(x::CurveDifference, t::Number) = value_at(x.a, t) - value_at(x.b, t)

function incremental_integral(x::AbstractCurve, t::Number, dt::Number)
    return integral(x, t + dt) - integral(x, t)
end

# The integral of the curve from 0 to x: trapezoids between the nodes, the rate flat before the first node and after the
# last one
function integral(curve::Curve, x::Number)
    T = curve.times
    r = curve.values
    iszero(x) && return zero(r[1] * x)
    ChainRulesCore.@ignore_derivatives check_domain(x >= 0, x, "The curve is integrated up to a nonnegative time")
    x <= T[1] && return r[1] * x
    idx_ = searchsortedlast(T, x; lt = <=)    # the last node before x
    @muladd out = r[1] * T[1] + sum((r[i] + r[i+1]) / 2 * (T[i+1] - T[i]) for i = 1:(idx_-1); init = zero(r[1] * x))
    @muladd return out + (r[idx_] + value_at(curve, x)) / 2 * (x - T[idx_])
end
integral(x::CurveDifference, t::Number) = integral(x.a, t) - integral(x.b, t)
integral(x::num1, t::num2) where {num1 <: Number, num2 <: Number} = x * t;

# Reverse rule of the integral of a curve (#80): Zygote's pullback of the element by element reads above cost ~0.5 ms
# per call (a dense one-hot array per getindex). The integral is r_1 min(x, T_1) + r_N (x - T_N)⁺ plus, for every
# segment i, the part of it before x: w (r_i + s w / 2), w = min(x, T_{i+1}) - T_i ≥ 0, s = (r_{i+1} - r_i) / h,
# h = T_{i+1} - T_i. Its partials, in closed form: in the values the trapezoid weights, in x the value of the curve, in
# the times through w and s. Only the segments before x contribute.
function integral_partials(curve::Curve{N}, x::Number) where {N}
    T, r = curve.times, curve.values
    FT = typeof(zero(eltype(T)) * zero(eltype(r)) * zero(x))
    dT, dr, dx = zeros(FT, N), zeros(FT, N), zero(FT)
    dr[1] = min(x, T[1])
    x < T[1] ? (dx += r[1]) : (dT[1] += r[1])
    if x > T[N]
        dr[N] += x - T[N]
        dx += r[N]
        dT[N] -= r[N]
    end
    for i = 1:(N-1)
        x <= T[i] && break
        h = T[i+1] - T[i]
        s = (r[i+1] - r[i]) / h
        inside = x < T[i+1]
        w = inside ? x - T[i] : h
        w2h = w^2 / (2 * h)
        dr[i] += w - w2h
        dr[i+1] += w2h
        @muladd v = r[i] + s * w  # the value of the curve at the end of the part
        inside && (dx += v)
        @muladd dT[i] += w2h * s - v
        @muladd dT[i+1] += ifelse(inside, zero(v), v) - w2h * s
    end
    return dT, dr, dx
end
function ChainRulesCore.rrule(::typeof(integral), curve::Curve{N}, x::Number) where {N}
    function integral_pullback(Δ)
        Δ_ = unthunk(Δ)
        dT, dr, dx = integral_partials(curve, x)
        return NoTangent(), Tangent{typeof(curve)}(; times = SVector{N}(Δ_ .* dT), values = SVector{N}(Δ_ .* dr)), Δ_ * dx
    end
    return integral(curve, x), integral_pullback
end
