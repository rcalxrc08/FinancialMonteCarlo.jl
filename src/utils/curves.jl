########### Curve Section

"""
    Curve(r::Vector, T::Vector)
    Curve(r::Vector, T)

Curve of instantaneous rates with values `r` at the sorted times `T` (or on a uniform grid from 0
to `T`), linearly interpolated in between and flat outside. Models only use its integrals over time
intervals, both for the drift and for discounting. Curves combine with `+`, `-`, `*` and `/`, with
numbers and with each other (on the union of their times).
"""
struct Curve{num1 <: Number, num2 <: Number}
    times::Vector{num1}
    values::Vector{num2}
    function Curve{num1, num2}(times::Vector{num1}, values::Vector{num2}) where {num1 <: Number, num2 <: Number}
        ChainRulesCore.@ignore_derivatives begin
            length(times) == length(values) || throw(ArgumentError("a curve needs one value per time"))
            isempty(times) && throw(ArgumentError("a curve needs at least one node"))
            issorted(times; lt = <=) || throw(ArgumentError("the times of a curve must be increasing"))
        end
        return new{num1, num2}(times, values)
    end
end
const CurveType = Curve
curve_from_nodes(times::AbstractVector, values::AbstractVector) = Curve{eltype(times), eltype(values)}(collect(times), collect(values))

function Curve(r_::AbstractVector{num1}, T::num2) where {num1 <: Number, num2 <: Number}
    N = length(r_)
    times = N == 1 ? [zero(T / 1)] : [(i - 1) * (T / (N - 1)) for i = 1:N]
    return curve_from_nodes(times, r_)
end
Curve(r_::AbstractVector{<:Number}, T::AbstractVector{<:Number}) = curve_from_nodes(T, r_)

"""
    ImpliedCurve(r::Vector, T::Vector)
    ImpliedCurve(r::Vector, T)

Piecewise linear curve of instantaneous rates, starting from zero at time zero, whose integral up to
each `T[i]` equals `r[i] * T[i]`: use it to turn zero rates (or dividend yields) quoted at maturities
`T` into a [`FinancialMonteCarlo.Curve`](@ref). With a scalar `T` the maturities are `T/N, 2T/N, ..., T`.
"""
function ImpliedCurve(r_::AbstractVector{num1}, T::AbstractVector{num2}) where {num1 <: Number, num2 <: Number}
    ChainRulesCore.@ignore_derivatives @assert length(r_) == length(T)
    ChainRulesCore.@ignore_derivatives @assert length(r_) >= 1
    # each rate is the one that makes the integral of the piecewise linear curve r_i T_i at T_i
    rates = [r_[1] * 2]
    for i = 2:length(r_)
        push!(rates, (r_[i] * T[i] - r_[i-1] * T[i-1]) * 2 / (T[i] - T[i-1]) - rates[end])
    end
    return curve_from_nodes(vcat(zero(num2), T), vcat(zero(eltype(rates)), rates))
end
ImpliedCurve(r_::AbstractVector{<:Number}, T::Number) = ImpliedCurve(r_, [i * (T / length(r_)) for i = 1:length(r_)])

extract_keys(curve::Curve) = curve.times
extract_values(curve::Curve) = curve.values
Base.keys(curve::Curve) = curve.times
Base.values(curve::Curve) = curve.values
Base.length(curve::Curve) = length(curve.times)
Base.map(f, curve::Curve) = curve_from_nodes(curve.times, map(f, curve.values))
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
    return r[i] + (t - T[i]) / (T[i+1] - T[i]) * (r[i+1] - r[i])
end

# Arithmetic: with a number on the values, between curves on the union of their times (each curve linear between its
# nodes and flat outside, #43, #44)
for op in (:+, :-, :*, :/)
    @eval Base.$op(x::Curve, y::Number) = map(value -> $op(value, y), x)
    @eval Base.$op(x::Number, y::Curve) = map(value -> $op(x, value), y)
    @eval function Base.$op(x::Curve, y::Curve)
        times = sort!(unique!(vcat(x.times, y.times)))
        return curve_from_nodes(times, [$op(value_at(x, t), value_at(y, t)) for t in times])
    end
end
Base.:-(curve::Curve) = map(-, curve)

function incremental_integral(x::Curve, t::Number, dt::Number)
    return integral(x, t + dt) - integral(x, t)
end

# The integral of the curve from 0 to x: trapezoids between the nodes, the rate flat before the first node and after the
# last one
function integral(curve::Curve, x::Number)
    T = curve.times
    r = curve.values
    iszero(x) && return zero(r[1] * x)
    ChainRulesCore.@ignore_derivatives @assert x >= 0
    x <= T[1] && return r[1] * x
    idx_ = searchsortedlast(T, x; lt = <=)    # the last node before x
    out = r[1] * T[1] + sum((r[i] + r[i+1]) / 2 * (T[i+1] - T[i]) for i = 1:(idx_-1); init = zero(r[1] * x))
    return out + (r[idx_] + value_at(curve, x)) / 2 * (x - T[idx_])
end
integral(x::num1, t::num2) where {num1 <: Number, num2 <: Number} = x * t;
