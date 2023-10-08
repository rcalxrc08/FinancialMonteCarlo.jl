"""
Class for Dispatching Forward Payoff

		forward=Forward(T::num) where {num <: Number}
	
Where:

		T	=	Time to maturity of the Forward.
"""
struct Forward{num <: Number} <: EuropeanPayoff{num}
    T::num
    function Forward(T::num) where {num <: Number}
        ChainRulesCore.@ignore_derivatives @assert T > 0 "Time to Maturity must be positive"
        return new{num}(T)
    end
end

export Forward;

function payoff(S::AbstractMatrix{num}, optionData::Forward, rfCurve::AbstractZeroRateCurve, mcBaseData::AbstractMonteCarloConfiguration, T1::num2 = maturity(optionData)) where {num <: Number, num2 <: Number}
    r = rfCurve.r
    T = optionData.T
    NStep = mcBaseData.Nstep
    index1 = round(Int, T / T1 * NStep) + 1
    @views ST = S[:, index1]

    return ST * exp(-integral(r, T))
end
