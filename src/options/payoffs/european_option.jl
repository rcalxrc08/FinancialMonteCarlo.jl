"""
    EuropeanOption(T, K, isCall = true)

European vanilla option paying ``\\max(S_T - K, 0)`` (call) or ``\\max(K - S_T, 0)`` (put) at maturity `T`.

# Arguments
- `T`: time to maturity, must be positive.
- `K`: strike, must be positive.
- `isCall`: `true` (default) for a call, `false` for a put.
"""
struct EuropeanOption{num1 <: Number, num2 <: Number, numtype <: Number} <: EuropeanPayoff{numtype}
    T::num1
    K::num2
    isCall::Bool
    function EuropeanOption(T::num1, K::num2, isCall::Bool = true) where {num1 <: Number, num2 <: Number}
        ChainRulesCore.@ignore_derivatives @assert T > 0 "Time to Maturity must be positive"
        ChainRulesCore.@ignore_derivatives @assert K > 0 "Strike Price must be positive"
        zero_typed = ChainRulesCore.@ignore_derivatives zero(num1) + zero(num2)
        return new{num1, num2, typeof(zero_typed)}(T, K, isCall)
    end
end

export EuropeanOption;

function payout_untyped(ST::numtype_, isCall::Bool, ::Type{<:EuropeanOption}) where {numtype_ <: Number}
    iscall = ChainRulesCore.@ignore_derivatives ifelse(isCall, Int8(1), Int8(-1))
    zero_typed = ChainRulesCore.@ignore_derivatives zero(ST)
    return max(iscall * (ST - 1), zero_typed)
end

payout_untyped(ST::Number, euPayoff::EuropeanOption) = payout_untyped(ST, euPayoff.isCall, EuropeanOption)

function payout_strike(ST::numtype_, K::Number, isCall::Bool, ::Type{<:EuropeanOption}) where {numtype_ <: Number}
    result = K * payout_untyped(ST / K, isCall, EuropeanOption)
    return result
end

payout(ST::Number, euPayoff::EuropeanOption) = payout_strike(ST, euPayoff.K, euPayoff.isCall, EuropeanOption)
