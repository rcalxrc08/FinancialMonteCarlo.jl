## Payoffs
abstract type AbstractPayoff{numtype <: Number} end

"""
    ControlVariates(variate, conf_variate, method = StandardMC())

Control variates method for [`MonteCarloConfiguration`](@ref). A pilot simulation with
configuration `conf_variate` estimates the optimal coefficient and the price of the control
`variate` (exact when `variate` is a [`Forward`](@ref)); the main simulation then uses `method`.
The control must have the same maturity as the priced payoff.
"""
struct ControlVariates{abstractPayoff <: AbstractPayoff{<:Number}, abstractMonteCarloConfiguration <: AbstractMonteCarloConfiguration, abstractMonteCarloMethod <: AbstractMonteCarloMethod} <: AbstractMonteCarloMethod
    variate::abstractPayoff
    conf_variate::abstractMonteCarloConfiguration
    curr_method::abstractMonteCarloMethod
    function ControlVariates(variate::T1, conf_variate::abstractMonteCarloConfiguration, curr_method::abstractMonteCarloMethod = StandardMC()) where {T1 <: AbstractPayoff{<:Number}, abstractMonteCarloConfiguration <: AbstractMonteCarloConfiguration, abstractMonteCarloMethod <: AbstractMonteCarloMethod}
        return new{T1, abstractMonteCarloConfiguration, abstractMonteCarloMethod}(variate, conf_variate, curr_method)
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
    payoff(S, payoff, rfCurve, mcConfig, T1 = maturity(payoff))

Discounted payout of `payoff` along each simulated path in `S`, i.e. the vector whose mean is the
price. `S` is the output of [`simulate`](@ref) up to time `T1`, which can be later than the
maturity of the payoff (as happens when several payoffs are priced on the same paths).

This is the generic interface: each payoff family has an engine implementing it on top of the
payoff specific `payout`, see the Extending section of the manual.
"""
function payoff(::AbstractMatrix, ::AbstractPayoff, ::AbstractZeroRateCurve, ::AbstractMonteCarloConfiguration, T1::Number = 0.0)
    error("Function used just for documentation")
end

function payoff(S::Array{Array{num, 2}, 1}, optionData::SingleNamePayoff, rfCurve::abstractZeroRateCurve, mcBaseData::AbstractMonteCarloConfiguration, T1::num2 = optionData.T) where {abstractZeroRateCurve <: AbstractZeroRateCurve, num <: Number, num2 <: Number}
    #switch to only in julia 1.4
    ChainRulesCore.@ignore_derivatives @assert length(S) == 1
    return payoff(S[1], optionData, rfCurve, mcBaseData, T1)
end
