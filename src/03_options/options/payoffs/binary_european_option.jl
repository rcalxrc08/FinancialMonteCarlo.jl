"""
    BinaryEuropeanOption(T, K, isCall = true)

European cash-or-nothing digital option paying 1 at maturity `T` if ``S_T > K`` (call) or ``S_T < K`` (put).

# Arguments
- `T`: time to maturity, must be positive.
- `K`: strike, must be positive.
- `isCall`: `true` (default) for a call, `false` for a put.
"""
@option_type BinaryEuropeanOption EuropeanPayoff T K

function payout_untyped(ST::numtype_, isCall::Bool, ::Type{<:BinaryEuropeanOption}) where {numtype_ <: Number}
    return ifelse((ST - 1) * call_put_sign(isCall) > 0, one(numtype_), zero(numtype_))
end

payout_untyped(ST::Number, euPayoff::BinaryEuropeanOption) = payout_untyped(ST, euPayoff.isCall, BinaryEuropeanOption)

payout_strike(ST::Number, K::Number, isCall::Bool, ::Type{<:BinaryEuropeanOption}) = payout_untyped(ST / K, isCall, BinaryEuropeanOption)

payout(ST::Number, euPayoff::BinaryEuropeanOption) = payout_strike(ST, euPayoff.K, euPayoff.isCall, BinaryEuropeanOption)
