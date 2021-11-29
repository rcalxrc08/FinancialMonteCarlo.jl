########### Curve Section
import Base.insert!
insert!(dic::Dict, k::Number, val::Number) = dic[k] = val;
return nothing;
using Dictionaries
#const FinMCDict=Dict
const FinMCDict = Dictionaries.Dictionary

#const Curve = FinMCDict{num1, num2} where {num1 <: Number, num2 <: Number}
const CurveType = FinMCDict{num1, num2} where {num1 <: Number, num2 <: Number}

function keys_(x::CurveType{num1, num2}) where {num1 <: Number, num2 <: Number}
    return sort(collect(keys(x)))
end

function values_(x::CurveType{num1, num2}) where {num1 <: Number, num2 <: Number}
    T_d = collect(keys(x))
    idx_ = sortperm(T_d)
    return collect(values(x))[idx_]
end

function Curve(r_::Array{num1}, T::num2) where {num1 <: Number, num2 <: Number}
    r = FinMCDict{num2, num1}()
    Nstep = length(r_) - 1
    dt = T / Nstep
    for i = 1:length(r_)
        insert!(r, (i - 1) * dt, r_[i])
    end
    return r
end
function Curve(r_::Array{num1}, T::Array{num2}) where {num1 <: Number, num2 <: Number}
    @assert length(r_) == length(T)
    @assert T == sort(T)
    r = FinMCDict{num2, num1}()
    Nstep = length(r_) - 1
    for i = 1:length(r_)
        insert!(r, T[i], r_[i])
    end
    return r
end
function ImpliedCurve(r_::Array{num1}, T::Array{num2}) where {num1 <: Number, num2 <: Number}
    @assert length(r_) == length(T)
    @assert length(r_) >= 1
    @assert T == sort(T)

    r = FinMCDict{num2, num1}()
    Nstep = length(r_) - 1
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
    dt_ = T / (N)
    t_ = collect(dt_:dt_:T)
    return ImpliedCurve(r_, t_)
end
function (x::CurveType)(t::Number, dt::Number)
    T = collect(keys_(x))
    r = collect(values_(x))
    return intgral_2(t + dt, T, r) - intgral_2(t, T, r)
end

using Interpolations
function intgral_2(x::num, T::Array{num1}, r::Array{num2}) where {num <: Number, num1 <: Number, num2 <: Number}
    @assert length(T) == length(r)
    if (x == 0.0)
        return 0.0
    end
    @assert x >= 0.0
    idx_ = findlast(y -> y < x, T)
    out = sum([(r[i] + r[i+1]) * 0.5 * (T[i+1] - T[i]) for i = 1:(idx_-1)])
    if x <= T[end]
        itp = LinearInterpolation([T[idx_], T[idx_+1]], [r[idx_], r[idx_+1]], extrapolation_bc = Flat())
        out = out + (r[idx_] + itp(x)) * 0.5 * (x - T[idx_])
    else
        #continuation
        out = out + (r[idx_] + r[idx_-1]) * 0.5 * (x - T[idx_])
    end
    return out
end

function integral(r::FinMCDict{num1, num2}, t::Number) where {num1 <: Number, num2 <: Number}
    T = collect(keys_(r))
    r = collect(values_(r))
    return intgral_2(t, T, r)
end
integral(x::num1, t::num2) where {num1 <: Number, num2 <: Number} = x * t;

### end curve section
