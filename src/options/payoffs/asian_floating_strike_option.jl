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

payout(S::AbstractArray{<:Number}, opt_::AsianFloatingStrikeOption) = vanilla_payout(S[end], mean(S), opt_)
