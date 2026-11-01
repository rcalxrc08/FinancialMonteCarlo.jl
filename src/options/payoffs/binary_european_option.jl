"""
    BinaryEuropeanOption(T, K, isCall = true)

European cash-or-nothing digital option paying 1 at maturity `T` if ``S_T > K`` (call) or ``S_T < K`` (put).

# Arguments
- `T`: time to maturity, must be positive.
- `K`: strike, must be positive.
- `isCall`: `true` (default) for a call, `false` for a put.
"""
struct BinaryEuropeanOption{num1 <: Number, num2 <: Number, numtype <: Number} <: EuropeanPayoff{numtype}
    T::num1
    K::num2
    isCall::Bool
    function BinaryEuropeanOption(T::num1, K::num2, isCall::Bool = true) where {num1 <: Number, num2 <: Number}
        ChainRulesCore.@ignore_derivatives @assert T > 0 "Time to Maturity must be positive"
        ChainRulesCore.@ignore_derivatives @assert K > 0 "Strike Price must be positive"
        zero_typed = ChainRulesCore.@ignore_derivatives zero(num1) + zero(num2)
        return new{num1, num2, typeof(zero_typed)}(T, K, isCall)
    end
end

export BinaryEuropeanOption;

function payout_untyped(ST::numtype_, isCall::Bool, ::Type{<:BinaryEuropeanOption}) where {numtype_ <: Number}
    iscall = ChainRulesCore.@ignore_derivatives ifelse(isCall, Int8(1), Int8(-1))
    res = ifelse((ST - 1) * iscall > 0, one(numtype_), zero(numtype_))
    return res
end

payout_untyped(ST::Number, euPayoff::BinaryEuropeanOption) = payout_untyped(ST, euPayoff.isCall, BinaryEuropeanOption)

payout_strike(ST::Number, K::Number, isCall::Bool, ::Type{<:BinaryEuropeanOption}) = payout_untyped(ST / K, isCall, BinaryEuropeanOption)

payout(ST::Number, euPayoff::BinaryEuropeanOption) = payout_strike(ST, euPayoff.K, euPayoff.isCall, BinaryEuropeanOption)
