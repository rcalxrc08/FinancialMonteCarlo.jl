crossed_below(S::AbstractArray, barrier::Number) = any(x -> x < barrier, S)
crossed_above(S::AbstractArray, barrier::Number) = any(x -> x > barrier, S)

payout(S::AbstractArray{<:Number}, payoff_::BarrierPayoff) = vanilla_payout(S[end], payoff_.K, payoff_) * is_alive(S, payoff_)

"""
    BarrierOptionDownOut(T, K, barrier, isCall = true)

Down-and-out barrier option: a European option with strike `K` and maturity `T` that is knocked out
if the underlying ever goes below `barrier`. The barrier is monitored on the simulation time grid.

# Arguments
- `T`: time to maturity, must be positive.
- `K`: strike, must be positive.
- `barrier`: lower barrier, must be positive.
- `isCall`: `true` (default) for a call, `false` for a put.
"""
@option_type BarrierOptionDownOut BarrierPayoff T K barrier

is_alive(S::AbstractArray, payoff_::BarrierOptionDownOut) = !crossed_below(S, payoff_.barrier)
