## Payoffs
abstract type AbstractPayoff{numtype <: Number} end
# A set of options of family T with the same maturity: not an AbstractPayoff, it has one value per strike
abstract type AbstractSmile{T} end

"""
    ControlVariates(variate, conf_variate, sampling = StandardMC())

Control variates sampling for [`MonteCarloConfiguration`](@ref). A pilot simulation with
configuration `conf_variate` estimates the optimal coefficient and the price of the control
`variate` (exact when `variate` is a [`Forward`](@ref)); the main simulation then uses `sampling`.
The control must have the same maturity as the priced payoff.
"""
struct ControlVariates{abstractPayoff <: AbstractPayoff{<:Number}, abstractMonteCarloConfiguration <: AbstractMonteCarloConfiguration, abstractSampling <: AbstractSampling} <: AbstractSampling
    variate::abstractPayoff
    conf_variate::abstractMonteCarloConfiguration
    sampling::abstractSampling
    function ControlVariates(variate::T1, conf_variate::abstractMonteCarloConfiguration, sampling::abstractSampling = StandardMC()) where {T1 <: AbstractPayoff{<:Number}, abstractMonteCarloConfiguration <: AbstractMonteCarloConfiguration, abstractSampling <: AbstractSampling}
        return new{T1, abstractMonteCarloConfiguration, abstractSampling}(variate, conf_variate, sampling)
    end
end
# The main simulation uses the sampling of the control variates
check_sampling(cv::ControlVariates, Nsim::Integer) = check_sampling(cv.sampling, Nsim)
check_sampling_mode(cv::ControlVariates, mode::AbstractMode) = check_sampling_mode(cv.sampling, mode)

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

maturity(x::AbstractArray) = maximum(maturity, x)

const POSITIVE_FIELD_MESSAGES = Dict(:T => "Time to Maturity must be positive", :K => "Strike Price must be positive", :barrier => "Barrier must be positive")

# Defines and exports an option with the given numeric fields, all positive, plus isCall. A docstring can precede it.
macro option_type(name, family, fields...)
    nums = [Symbol(:num, i) for i in eachindex(fields)]
    typed_fields = [:($field::$num) for (field, num) in zip(fields, nums)]
    bounds = [:($num <: Number) for num in nums]
    checks = [:(ChainRulesCore.@ignore_derivatives @assert $field > 0 $(POSITIVE_FIELD_MESSAGES[field])) for field in fields]
    zeros_ = [:(zero($num)) for num in nums]
    return esc(quote
        Base.@__doc__ struct $name{$(bounds...), numtype <: Number} <: $family{numtype}
            $(typed_fields...)
            isCall::Bool
            function $name($(typed_fields...), isCall::Bool = true) where {$(bounds...)}
                $(checks...)
                zero_typed = ChainRulesCore.@ignore_derivatives +($(zeros_...))
                return new{$(nums...), typeof(zero_typed)}($(fields...), isCall)
            end
        end
        export $name
    end)
end

call_put_sign(isCall::Bool) = ChainRulesCore.@ignore_derivatives ifelse(isCall, Int8(1), Int8(-1))

# The value of a vanilla option, of the type of all the inputs of the payoff, whichever enter the formula
function vanilla_payout(x::Number, K::Number, payoff_::AbstractPayoff{numtype}) where {numtype <: Number}
    zero_typed = ChainRulesCore.@ignore_derivatives zero(x) * zero(numtype)
    return max(call_put_sign(payoff_.isCall) * (x - K), zero_typed)
end

# The column of the paths, simulated up to T1, of the maturity T
maturity_index(T::Number, T1::Number, mcBaseData::AbstractMonteCarloConfiguration) = round(Int, T / T1 * mcBaseData.Nstep) + 1

"""
    payoff(S, payoff, rfCurve, mcConfig, T1 = maturity(payoff))

Discounted payout of `payoff` along each simulated path in `S`, i.e. the vector whose mean is the
price. `S` is the output of [`simulate`](@ref) up to time `T1`, which can be later than the
maturity of the payoff (as happens when several payoffs are priced on the same paths).

This is the generic interface: each payoff family has an engine implementing it on top of the
payoff specific `payout`, see the Extending section of the manual.
"""
function payoff(::AbstractMatrix, ::AbstractPayoff, ::AbstractZeroRateCurve, ::AbstractMonteCarloConfiguration, T1::Number = 0)
    error("Function used just for documentation")
end

function payoff(S::Array{Array{num, 2}, 1}, optionData::Union{SingleNamePayoff, AbstractSmile}, rfCurve::abstractZeroRateCurve, mcBaseData::AbstractMonteCarloConfiguration, T1::num2 = optionData.T) where {abstractZeroRateCurve <: AbstractZeroRateCurve, num <: Number, num2 <: Number}
    #switch to only in julia 1.4
    ChainRulesCore.@ignore_derivatives @assert length(S) == 1
    return payoff(S[1], optionData, rfCurve, mcBaseData, T1)
end
