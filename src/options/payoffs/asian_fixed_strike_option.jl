"""
    AsianFixedStrikeOption(T, K, isCall = true)

Arithmetic Asian option with fixed strike, paying ``\\max(A_T - K, 0)`` (call) or ``\\max(K - A_T, 0)``
(put) at maturity `T`, where ``A_T`` is the average of the underlying on the simulation time grid.

# Arguments
- `T`: time to maturity, must be positive.
- `K`: strike, must be positive.
- `isCall`: `true` (default) for a call, `false` for a put.
"""
struct AsianFixedStrikeOption{num1 <: Number, num2 <: Number, numtype <: Number} <: AsianPayoff{numtype}
    T::num1
    K::num2
    isCall::Bool
    function AsianFixedStrikeOption(T::num1, K::num2, isCall::Bool = true) where {num1 <: Number, num2 <: Number}
        ChainRulesCore.@ignore_derivatives @assert T > 0 "Time to Maturity must be positive"
        ChainRulesCore.@ignore_derivatives @assert K > 0 "Strike Price must be positive"
        zero_typed = ChainRulesCore.@ignore_derivatives zero(num1) + zero(num2)
        return new{num1, num2, typeof(zero_typed)}(T, K, isCall)
    end
end

export AsianFixedStrikeOption;

function payout(S::abstractArray, opt_::AsianFixedStrikeOption) where {abstractArray <: AbstractArray{num_}} where {num_ <: Number}
    iscall = ChainRulesCore.@ignore_derivatives ifelse(opt_.isCall, Int8(1), Int8(-1))
    zero_typed = ChainRulesCore.@ignore_derivatives zero(eltype(S)) * opt_.K
    return max(iscall * (mean(S) - opt_.K), zero_typed)
end
