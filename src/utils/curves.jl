########### Curve Section
using Dictionaries

const DictTypeInternal = Dictionaries.Dictionary

const CurveType = DictTypeInternal{num1, num2} where {num1 <: Number, num2 <: Number}

function extract_keys(x::CurveType{num1, num2}) where {num1 <: Number, num2 <: Number}
    return sort(collect(keys(x)))
end

function extract_values(x::CurveType{num1, num2}) where {num1 <: Number, num2 <: Number}
    T_d = collect(keys(x))
    idx_ = sortperm(T_d)
    return collect(values(x))[idx_]
end

function Curve(r_::Array{num1}, T::num2) where {num1 <: Number, num2 <: Number}
    r = DictTypeInternal{num2, num1}()
    Nstep = length(r_) - 1
    dt = T / Nstep
    for i = 1:length(r_)
        insert!(r, (i - 1) * dt, r_[i])
    end
    return r
end
function Curve(r_::Array{num1}, T::Array{num2}) where {num1 <: Number, num2 <: Number}
    ChainRulesCore.@ignore_derivatives @assert length(r_) == length(T)
    ChainRulesCore.@ignore_derivatives @assert T == sort(T)
    r = DictTypeInternal{num2, num1}()
    for i = 1:length(r_)
        insert!(r, T[i], r_[i])
    end
    return r
end
function ImpliedCurve(r_::Array{num1}, T::Array{num2}) where {num1 <: Number, num2 <: Number}
    ChainRulesCore.@ignore_derivatives @assert length(r_) == length(T)
    ChainRulesCore.@ignore_derivatives @assert length(r_) >= 1
    ChainRulesCore.@ignore_derivatives @assert T == sort(T)

    r = DictTypeInternal{num2, num1}()
    insert!(r, 0.0, 0.0)
    insert!(r, T[1], r_[1] * 2.0)
    prec_r = r_[1] * 2.0
    for i = 2:length(r_)
        tmp_r = (r_[i] * T[i] - r_[i-1] * T[i-1]) * 2 / (T[i] - T[i-1]) - prec_r
        insert!(r, T[i], tmp_r)
        prec_r = tmp_r
    end
    return r
end
function ImpliedCurve(r_::Array{num1}, T::num2) where {num1 <: Number, num2 <: Number}
    N = length(r_)
    dt_ = T / N
    t_ = collect(dt_:dt_:T)
    return ImpliedCurve(r_, t_)
end

function incremental_integral(x::CurveType{num1, num2}, t::Number, dt::Number) where {num1 <: Number, num2 <: Number}
    T = extract_keys(x)
    r = extract_values(x)
    return internal_definite_integral(t + dt, T, r) - internal_definite_integral(t, T, r)
end

using Interpolations
function internal_definite_integral(x::num, T::Array{num1}, r::Array{num2}) where {num <: Number, num1 <: Number, num2 <: Number}
    ChainRulesCore.@ignore_derivatives @assert length(T) == length(r)
    if x == 0.0
        return 0.0
    end
    ChainRulesCore.@ignore_derivatives @assert x >= 0.0
    idx_ = findlast(y -> y < x, T)
    out = sum([(r[i] + r[i+1]) * 0.5 * (T[i+1] - T[i]) for i = 1:(idx_-1)])
    if x <= T[end]
        itp = linear_interpolation([T[idx_], T[idx_+1]], [r[idx_], r[idx_+1]], extrapolation_bc = Flat())
        out = out + (r[idx_] + itp(x)) * 0.5 * (x - T[idx_])
    else
        #continuation
        out = out + (r[idx_] + r[idx_-1]) * 0.5 * (x - T[idx_])
    end
    return out
end

function integral(r::DictTypeInternal{num1, num2}, t::Number) where {num1 <: Number, num2 <: Number}
    T = extract_keys(r)
    values = extract_values(r)
    return internal_definite_integral(t, T, values)
end
integral(x::num1, t::num2) where {num1 <: Number, num2 <: Number} = x * t;

### end curve section
