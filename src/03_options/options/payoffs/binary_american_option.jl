"""
    BinaryAmericanOption(T, K, isCall = true)

American digital option, exercisable at any time step up to maturity `T`, paying 1 when exercised
in the money (``S_t > K`` for a call, ``S_t < K`` for a put). Priced with the Longstaff-Schwartz regression.

# Arguments
- `T`: time to maturity, must be positive.
- `K`: strike, must be positive.
- `isCall`: `true` (default) for a call, `false` for a put.
"""
@option_type BinaryAmericanOption AmericanPayoff T K

function payout(Sti::numtype_, amPayoff::BinaryAmericanOption) where {numtype_ <: Number}
    return ifelse((Sti - amPayoff.K) * call_put_sign(amPayoff.isCall) > 0, one(numtype_), zero(numtype_))
end
