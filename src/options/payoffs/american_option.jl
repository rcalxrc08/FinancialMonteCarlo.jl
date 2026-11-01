"""
    AmericanOption(T, K, isCall = true)

American vanilla option, exercisable at any time step up to maturity `T` with exercise value
``\\max(S_t - K, 0)`` (call) or ``\\max(K - S_t, 0)`` (put). Priced with the Longstaff-Schwartz regression.

# Arguments
- `T`: time to maturity, must be positive.
- `K`: strike, must be positive.
- `isCall`: `true` (default) for a call, `false` for a put.
"""
@option_type AmericanOption AmericanPayoff T K

payout(Sti::Number, amPayoff::AmericanOption) = vanilla_payout(Sti, amPayoff.K, amPayoff)
