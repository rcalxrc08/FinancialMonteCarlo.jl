## Payoffs
abstract type AbstractPayoff end

struct ControlVariates{abstractPayoff <: AbstractPayoff} <: AbstractMonteCarloMethod
	variate::AbstractPayoff
	conf_variate::MonteCarloConfiguration{<: Integer , <: Integer , <: AbstractMonteCarloMethod ,  <: BaseMode}
	function ControlVariates(variate::T1,conf_variate::MonteCarloConfiguration) where { T1 <: AbstractPayoff}
		return new{T1}(variate,conf_variate)
    end

end

abstract type SingleNamePayoff<:AbstractPayoff end
abstract type EuropeanPayoff<:SingleNamePayoff end
abstract type AmericanPayoff<:SingleNamePayoff end
abstract type PathDependentPayoff<:SingleNamePayoff end
abstract type BarrierPayoff<:PathDependentPayoff end
abstract type AsianPayoff<:PathDependentPayoff end

#No multiple inheritance, sigh
abstract type BasketPayoff<:AbstractPayoff end
abstract type EuropeanBasketPayoff<:BasketPayoff end
abstract type AmericanBasketPayoff<:BasketPayoff end
abstract type BarrierBasketPayoff<:BasketPayoff end
abstract type AsianBasketPayoff<:BasketPayoff end

####### Payoffs definition

### Spot
include("payoffs/spot.jl")

### European Payoffs
include("payoffs/general_european_engine.jl")
include("payoffs/forward.jl")
include("payoffs/european_option.jl")
include("payoffs/binary_european_option.jl")

### Barrier Payoffs
include("payoffs/barrier_do_option.jl")
include("payoffs/barrier_di_option.jl")
include("payoffs/barrier_uo_option.jl")
include("payoffs/barrier_ui_option.jl")
include("payoffs/double_barrier_option.jl")

### American Payoffs
include("payoffs/general_american_engine.jl")
include("payoffs/american_option.jl")
include("payoffs/binary_american_option.jl")

### Path Dependent Payoffs
include("payoffs/general_path_dependent_engine.jl")

### Asian Payoffs
include("payoffs/asian_fixed_strike_option.jl")
include("payoffs/asian_floating_strike_option.jl")

### Basket Payoffs
include("payoffs/basket/n_european_option.jl")

### Operations and Strategies
include("payoffs/operations.jl")
include("payoffs/operations2.jl")

##Support

function maturity(x::abstractPayoff) where { abstractPayoff <: AbstractPayoff}
	return x.T;
end

function maturity(x::Spot)
	return 0.0;
end


"""
General Interface for Payoff computation from MonteCarlo paths

		Payoff=payoff(S,optionData,rfCurve,T1=optionData.T)
	
Where:\n
		S           = Paths of the Underlying.
		optionData  = Datas of the Option.
		rfCurve  = Datas of the Risk Free Curve.
		T1 [Optional, default to optionData.T]=Final Time of Spot Simulation (default equals Time to Maturity of the option)

		Payoff      = payoff of the Option.

"""
function payoff(S::AbstractMatrix,optionData::AbstractPayoff,rfCurve::AbstractZeroRateCurve,T1::Number=optionData.T)
	error("Function used just for documentation")
end

function payoff(S::Array{Array{num,2},1},optionData::SingleNamePayoff,rfCurve::abstractZeroRateCurve,T1::num2=optionData.T) where{ abstractZeroRateCurve <: AbstractZeroRateCurve,num <: Number,num2 <: Number}
	@assert length(S)==1
	return payoff(S[1],optionData,rfCurve,T1);
end

