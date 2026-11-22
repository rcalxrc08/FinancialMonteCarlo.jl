"""
    DoubleBarrierOption(T, K, D, U, isCall = true)

Double knock-out barrier option: a European option with strike `K` and maturity `T` that is knocked
out if the underlying ever goes below `D` or above `U`. The barrier is monitored on the simulation time grid.

# Arguments
- `T`: time to maturity, must be positive.
- `K`: strike, must be positive.
- `D`: lower barrier, must be positive.
- `U`: upper barrier, must be positive.
- `isCall`: `true` (default) for a call, `false` for a put.
"""
struct DoubleBarrierOption{num1 <: Number, num2 <: Number, num3 <: Number, num4 <: Number, numtype <: Number} <: BarrierPayoff{numtype}
    T::num1
    K::num2
    D::num3
    U::num4
    isCall::Bool
    function DoubleBarrierOption(T::num1, K::num2, D::num3, U::num4, isCall::Bool = true) where {num1 <: Number, num2 <: Number, num3 <: Number, num4 <: Number}
        check_double_barrier(T, K, D, U)
        zero_typed = zero_of_types(num1, num2, num3, num4)
        return new{num1, num2, num3, num4, typeof(zero_typed)}(T, K, D, U, isCall)
    end
end

function check_double_barrier(T, K, D, U)
    check_maturity(T)
    check_domain(K > 0, K, "Strike Price must be positive")
    check_domain(D > 0, D, "Low Barrier must be positive")
    check_domain(U > D, U, "High Barrier must be positive and greater then Low barrier")
end
ChainRulesCore.@non_differentiable check_double_barrier(::Any, ::Any, ::Any, ::Any)

export DoubleBarrierOption;

is_alive(S, payoff_::DoubleBarrierOption) = !crossed_below(S, payoff_.D) && !crossed_above(S, payoff_.U)
