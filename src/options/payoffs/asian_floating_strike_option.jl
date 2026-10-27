"""
    AsianFloatingStrikeOption(T, isCall = true)

Arithmetic Asian option with floating strike, paying ``\\max(S_T - A_T, 0)`` (call) or ``\\max(A_T - S_T, 0)``
(put) at maturity `T`, where ``A_T`` is the average of the underlying on the simulation time grid.

# Arguments
- `T`: time to maturity, must be positive.
- `isCall`: `true` (default) for a call, `false` for a put.
"""
struct AsianFloatingStrikeOption{num <: Number} <: AsianPayoff{num}
    T::num
    isCall::Bool
    function AsianFloatingStrikeOption(T::num, isCall::Bool = true) where {num <: Number}
        ChainRulesCore.@ignore_derivatives @assert T > 0 "Time to Maturity must be positive"
        return new{num}(T, isCall)
    end
end

export AsianFloatingStrikeOption;

function payout(S::abstractArray, opt_::AsianFloatingStrikeOption) where {abstractArray <: AbstractArray{num_}} where {num_ <: Number}
    iscall = ChainRulesCore.@ignore_derivatives ifelse(opt_.isCall, Int8(1), Int8(-1))
    zero_typed = ChainRulesCore.@ignore_derivatives zero(eltype(S))
    return max(iscall * (S[end] - mean(S)), zero_typed)
end
