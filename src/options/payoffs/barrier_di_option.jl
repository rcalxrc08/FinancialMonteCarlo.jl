"""
    BarrierOptionDownIn(T, K, barrier, isCall = true)

Down-and-in barrier option: a European option with strike `K` and maturity `T` that pays only
if the underlying goes below `barrier` at least once. The barrier is monitored on the simulation time grid.

# Arguments
- `T`: time to maturity, must be positive.
- `K`: strike, must be positive.
- `barrier`: lower barrier, must be positive.
- `isCall`: `true` (default) for a call, `false` for a put.
"""
@option_type BarrierOptionDownIn BarrierPayoff T K barrier

is_alive(S::AbstractArray, payoff_::BarrierOptionDownIn) = crossed_below(S, payoff_.barrier)
