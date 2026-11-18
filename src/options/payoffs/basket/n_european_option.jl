"""
    EuropeanOptionND(T, K, isCall = true)

European option on the sum of the underlyings of a multivariate model, paying
``\\max(\\sum_i S^i_T - K, 0)`` (call) or ``\\max(K - \\sum_i S^i_T, 0)`` (put) at maturity `T`.

# Arguments
- `T`: time to maturity, must be positive.
- `K`: strike, must be positive.
- `isCall`: `true` (default) for a call, `false` for a put.
"""
struct EuropeanOptionND{num1 <: Number, num2 <: Number, numtype <: Number} <: EuropeanBasketPayoff{numtype}
    T::num1
    K::num2
    isCall::Bool
    function EuropeanOptionND(T::num1, K::num2, isCall::Bool = true) where {num1 <: Number, num2 <: Number}
        ChainRulesCore.@ignore_derivatives check_domain(T > 0, T, "Time to Maturity must be positive")
        ChainRulesCore.@ignore_derivatives check_domain(K > 0, K, "Strike Price must be positive")
        zero_typed = ChainRulesCore.@ignore_derivatives zero(num1) + zero(num2)
        return new{num1, num2, typeof(zero_typed)}(T, K, isCall)
    end
end

export EuropeanOptionND;

function payoff(S::Array{abstractMatrix}, euPayoff::EuropeanOptionND, rfCurve::abstractZeroRateCurve, mcBaseData::AbstractMonteCarloConfiguration, T1::num2 = maturity(euPayoff)) where {abstractMatrix <: AbstractMatrix, num2 <: Number, abstractZeroRateCurve <: AbstractZeroRateCurve}
    r = rfCurve.r
    T = euPayoff.T
    iscall = ChainRulesCore.@ignore_derivatives ifelse(euPayoff.isCall, Int8(1), Int8(-1))
    NStep = mcBaseData.Nstep
    index1 = maturity_index(T, T1, mcBaseData)
    K = euPayoff.K
    #ST_all=[ sum(x_i[j,index1] for x_i in S) for j in 1:Nsim]
    # ST_all = S[1][:, index1]
    # for i = 2:length(S)
    #     ST_all += S[i][:, index1]
    # end
    ST_all = sum(view(S[i], :, index1) for i = 1:length(S))
    # a zero of the type of the payout: zero(K) alone is a Float64 next to paths of hyper duals, and max(::Hyper, ::Float64)
    # returns either argument unconverted, so the vector would get an abstract element type (#53)
    zero_typed = ChainRulesCore.@ignore_derivatives zero(eltype(ST_all)) - zero(K)
    payoff2 = max.(iscall .* (ST_all .- K), zero_typed)

    return payoff2 * exp(-integral(r, T))
end
