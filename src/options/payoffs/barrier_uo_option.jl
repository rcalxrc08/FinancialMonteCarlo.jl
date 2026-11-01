"""
    BarrierOptionUpOut(T, K, barrier, isCall = true)

Up-and-out barrier option: a European option with strike `K` and maturity `T` that is knocked out
if the underlying ever goes above `barrier`. The barrier is monitored on the simulation time grid.

# Arguments
- `T`: time to maturity, must be positive.
- `K`: strike, must be positive.
- `barrier`: upper barrier, must be positive.
- `isCall`: `true` (default) for a call, `false` for a put.
"""
@option_type BarrierOptionUpOut BarrierPayoff T K barrier

is_alive(S::AbstractArray, payoff_::BarrierOptionUpOut) = !crossed_above(S, payoff_.barrier)
