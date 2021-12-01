"""
Struct for European Option

		euOption=EuropeanOptionND(T::num1,K::num2,isCall::Bool=true) where {num1 <: Number,num2 <: Number}
	
Where:\n
		T	=	Time to maturity of the Option.
		K	=	Strike Price of the Option.
		isCall  = true for CALL, false for PUT.
"""
mutable struct EuropeanOptionND{num1 <: Number, num2 <: Number, numtype <: Number} <: EuropeanBasketPayoff{numtype}
    T::num1
    K::num2
    isCall::Bool
    function EuropeanOptionND(T::num1, K::num2, isCall::Bool = true) where {num1 <: Number, num2 <: Number}
        @assert T > 0 "Time to Maturity must be positive"
        @assert K > 0 "Strike Price must be positive"
        zero_typed = zero(num1) + zero(num2)
        return new{num1, num2, typeof(zero_typed)}(T, K, isCall)
    end
end

export EuropeanOptionND;

function payoff(S::Array{abstractMatrix}, euPayoff::EuropeanOptionND, rfCurve::abstractZeroRateCurve, T1::num2 = maturity(euPayoff)) where {abstractMatrix <: AbstractMatrix, num2 <: Number, abstractZeroRateCurve <: AbstractZeroRateCurve}
    r = rfCurve.r
    T = euPayoff.T
    iscall = euPayoff.isCall ? 1 : -1
    (Nsim, NStep) = size(S[1])
    NStep -= 1
    index1 = round(Int, T / T1 * NStep) + 1
    K = euPayoff.K
    #ST_all=[ sum(x_i[j,index1] for x_i in S) for j in 1:Nsim]
    # ST_all = S[1][:, index1]
    # for i = 2:length(S)
    #     ST_all += S[i][:, index1]
    # end
    ST_all = sum(S[i][:, index1] for i = 1:length(S))
    payoff2 = max.(iscall * (ST_all .- K), 0.0)

    return payoff2 * exp(-integral(r, T))
end
