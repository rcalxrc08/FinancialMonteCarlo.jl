"""
    AsianFixedStrikeOption(T, K, isCall = true)

Arithmetic Asian option with fixed strike, paying ``\\max(A_T - K, 0)`` (call) or ``\\max(K - A_T, 0)``
(put) at maturity `T`, where ``A_T`` is the average of the underlying on the simulation time grid.

# Arguments
- `T`: time to maturity, must be positive.
- `K`: strike, must be positive.
- `isCall`: `true` (default) for a call, `false` for a put.
"""
@option_type AsianFixedStrikeOption AsianPayoff T K

payout(S::AbstractArray{<:Number}, opt_::AsianFixedStrikeOption) = vanilla_payout(mean(S), opt_.K, opt_)
# The payout from the final value and the mean of the path (payoff_by_rows)
asian_payout(ST::Number, A::Number, opt_::AsianFixedStrikeOption) = vanilla_payout(A, opt_.K, opt_)
