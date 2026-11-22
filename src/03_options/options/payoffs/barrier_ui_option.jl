"""
    BarrierOptionUpIn(T, K, barrier, isCall = true)

Up-and-in barrier option: a European option with strike `K` and maturity `T` that pays only
if the underlying goes above `barrier` at least once. The barrier is monitored on the simulation time grid.

# Arguments
- `T`: time to maturity, must be positive.
- `K`: strike, must be positive.
- `barrier`: upper barrier, must be positive.
- `isCall`: `true` (default) for a call, `false` for a put.
"""
@option_type BarrierOptionUpIn BarrierPayoff T K barrier

is_alive(S, payoff_::BarrierOptionUpIn) = crossed_above(S, payoff_.barrier)
