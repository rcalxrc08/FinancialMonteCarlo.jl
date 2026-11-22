"""
    BermudanOption(T, K, isCall = true)

Bermudan vanilla option, exercisable at the dates in `T` (the last one being the maturity), with
exercise value ``\\max(S_t - K, 0)`` (call) or ``\\max(K - S_t, 0)`` (put). Priced with the
Longstaff-Schwartz regression. Exercise dates should lie on the simulation time grid.

# Arguments
- `T`: sorted vector of positive exercise dates.
- `K`: strike, must be positive.
- `isCall`: `true` (default) for a call, `false` for a put.
"""
struct BermudanOption{num1 <: Number, num2 <: Number, numtype <: Number} <: BermudanPayoff{numtype}
    T::num1
    T_ex::Vector{num1}
    K::num2
    isCall::Bool
    function BermudanOption(T::AbstractVector{num1}, K::num2, isCall::Bool = true) where {num1 <: Number, num2 <: Number}
        check_bermudan(T, K)
        zero_typed = zero_of_types(num1, num2)
        return new{num1, num2, typeof(zero_typed)}(T[end], T, K, isCall)
    end
end

function check_bermudan(T, K)
    check_domain(minimum(T) > 0, T, "Times to Maturity must be positive")
    check_domain(K > 0, K, "Strike Price must be positive")
    check_argument(issorted(T), "Times to Maturity must be sorted")
end
ChainRulesCore.@non_differentiable check_bermudan(::Any, ::Any)

export BermudanOption;

payout(Sti::Number, bmPayoff::BermudanOption) = vanilla_payout(Sti, bmPayoff.K, bmPayoff)
