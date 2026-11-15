"""
Struct for Average Forward

		forward=ForwardND(T::num1) where {num1 <: Number}
	
Where:\n
		T	=	Time to maturity of the Option.
"""
struct ForwardND{num1 <: Number} <: EuropeanBasketPayoff{num1}
    T::num1
    function ForwardND(T::num1) where {num1 <: Number}
        ChainRulesCore.@ignore_derivatives @assert T > 0 "Time to Maturity must be positive"
        return new{num1}(T)
    end
end

export ForwardND;

function payoff(S::Array{abstractMatrix}, euPayoff::ForwardND, rfCurve::abstractZeroRateCurve, T1::num2 = maturity(euPayoff)) where {abstractMatrix <: AbstractMatrix, num2 <: Number, abstractZeroRateCurve <: AbstractZeroRateCurve}
    r = rfCurve.r
    T = euPayoff.T
    (Nsim, NStep) = size(S[1])
    NStep -= 1
    index1 = round(Int, T / T1 * NStep) + 1
    #ST_all=[ sum(x_i[j,index1] for x_i in S) for j in 1:Nsim]
    ST_all = S[1][:, index1]
    for i = 2:length(S)
        ST_all += S[i][:, index1]
    end

    return ST_all * exp(-integral(r, T)) / length(S)
end
