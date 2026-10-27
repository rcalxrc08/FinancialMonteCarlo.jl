"""
    BermudanOption(T, K, isCall = true)

Bermudan vanilla option, exercisable at the dates in `T` (the last one being the maturity), with
exercise value ``\\max(S_t - K, 0)`` (call) or ``\\max(K - S_t, 0)`` (put). Priced with the
Longstaff-Schwartz regression. Exercise dates should lie on the simulation time grid.

# Arguments
- `T`: sorted vector of positive exercise dates.
- `K`: strike, must be positive.
- `isCall`: `true` (default) for a call, `false` for a put.
"""
struct BermudanOption{num1 <: Number, num2 <: Number, numtype <: Number} <: BermudanPayoff{numtype}
    T::num1
    T_ex::Array{num1}
    K::num2
    isCall::Bool
    function BermudanOption(T::AbstractArray{num1}, K::num2, isCall::Bool = true) where {num1 <: Number, num2 <: Number}
        ChainRulesCore.@ignore_derivatives @assert minimum(T) > 0 "Times to Maturity must be positive"
        ChainRulesCore.@ignore_derivatives @assert K > 0 "Strike Price must be positive"
        ChainRulesCore.@ignore_derivatives @assert issorted(T) "Times to Maturity must be sorted"
        zero_typed = ChainRulesCore.@ignore_derivatives zero(num1) + zero(num2)
        return new{num1, num2, typeof(zero_typed)}(T[end], T, K, isCall)
    end
end

export BermudanOption;

function payout(Sti::numtype_, bmPayoff::BermudanOption) where {numtype_ <: Number}
    iscall = ChainRulesCore.@ignore_derivatives ifelse(bmPayoff.isCall, Int8(1), Int8(-1))
    return max((Sti - bmPayoff.K) * iscall, zero(numtype_))
end
