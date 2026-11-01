"""
    EuropeanOption(T, K, isCall = true)

European vanilla option paying ``\\max(S_T - K, 0)`` (call) or ``\\max(K - S_T, 0)`` (put) at maturity `T`.

# Arguments
- `T`: time to maturity, must be positive.
- `K`: strike, must be positive.
- `isCall`: `true` (default) for a call, `false` for a put.
"""
@option_type EuropeanOption EuropeanPayoff T K

function payout_untyped(ST::numtype_, isCall::Bool, ::Type{<:EuropeanOption}) where {numtype_ <: Number}
    zero_typed = ChainRulesCore.@ignore_derivatives zero(ST)
    return max(call_put_sign(isCall) * (ST - 1), zero_typed)
end

payout_untyped(ST::Number, euPayoff::EuropeanOption) = payout_untyped(ST, euPayoff.isCall, EuropeanOption)

function payout_strike(ST::numtype_, K::Number, isCall::Bool, ::Type{<:EuropeanOption}) where {numtype_ <: Number}
    result = K * payout_untyped(ST / K, isCall, EuropeanOption)
    return result
end

payout(ST::Number, euPayoff::EuropeanOption) = payout_strike(ST, euPayoff.K, euPayoff.isCall, EuropeanOption)
