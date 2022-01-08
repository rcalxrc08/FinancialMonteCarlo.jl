## Payoffs
abstract type AbstractPayoff{numtype <: Number} end

struct ControlVariates{abstractPayoff <: AbstractPayoff{<:Number}} <: AbstractMonteCarloMethod
    variate::AbstractPayoff
    conf_variate::MonteCarloConfiguration{<:Integer, <:Integer, <:AbstractMonteCarloMethod, <:BaseMode}
    function ControlVariates(variate::T1, conf_variate::MonteCarloConfiguration) where {T1 <: AbstractPayoff{<:Number}}
        return new{T1}(variate, conf_variate)
    end
end

abstract type SingleNamePayoff{numtype <: Number} <: AbstractPayoff{numtype} end
abstract type EuropeanPayoff{numtype <: Number} <: SingleNamePayoff{numtype} end
abstract type AmericanPayoff{numtype <: Number} <: SingleNamePayoff{numtype} end
abstract type BermudanPayoff{numtype <: Number} <: SingleNamePayoff{numtype} end
abstract type PathDependentPayoff{numtype <: Number} <: SingleNamePayoff{numtype} end
abstract type BarrierPayoff{numtype <: Number} <: PathDependentPayoff{numtype} end
abstract type AsianPayoff{numtype <: Number} <: PathDependentPayoff{numtype} end

#No multiple inheritance, sigh
abstract type BasketPayoff{numtype <: Number} <: AbstractPayoff{numtype} end
abstract type EuropeanBasketPayoff{numtype <: Number} <: BasketPayoff{numtype} end
abstract type AmericanBasketPayoff{numtype <: Number} <: BasketPayoff{numtype} end
abstract type BarrierBasketPayoff{numtype <: Number} <: BasketPayoff{numtype} end
abstract type AsianBasketPayoff{numtype <: Number} <: BasketPayoff{numtype} end

function maturity(x::abstractPayoff) where {abstractPayoff <: AbstractPayoff{<:Number}}
    return x.T
end

"""
General Interface for Payoff computation from MonteCarlo paths

		Payoff=payoff(S,optionData,rfCurve,T1=optionData.T)
	
Where:

		S           = Paths of the Underlying.
		optionData  = Datas of the Option.
		rfCurve     = Datas of the Risk Free Curve.
		T1          =Final Time of Spot Simulation (default equals Time to Maturity of the option)[Optional, default to optionData.T]

		Payoff      = payoff of the Option.

"""
function payoff(::AbstractMatrix, ::AbstractPayoff, ::AbstractZeroRateCurve, ::AbstractMonteCarloConfiguration, T1::Number = 0.0)
    error("Function used just for documentation")
end

function payoff(S::Array{Array{num, 2}, 1}, optionData::SingleNamePayoff, rfCurve::abstractZeroRateCurve, mcBaseData::AbstractMonteCarloConfiguration, T1::num2 = optionData.T) where {abstractZeroRateCurve <: AbstractZeroRateCurve, num <: Number, num2 <: Number}
    #switch to only in julia 1.4
    @assert length(S) == 1
    return payoff(S[1], optionData, rfCurve, mcBaseData, T1)
end
