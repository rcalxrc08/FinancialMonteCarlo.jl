## Payoffs
abstract type AbstractPayoff end

abstract type EuropeanPayoff<:AbstractPayoff end
abstract type AmericanPayoff<:AbstractPayoff end
abstract type BarrierPayoff<:AbstractPayoff end
abstract type AsianPayoff<:AbstractPayoff end

####### Payoffs definition

### European Payoffs
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
include("payoffs/general_american_option.jl")
include("payoffs/american_option.jl")
include("payoffs/binary_american_option.jl")

### Asian Payoffs
include("payoffs/asian_fixed_strike_option.jl")
include("payoffs/asian_floating_strike_option.jl")


"""
General Interface for Payoff computation from MonteCarlo paths

		Payoff=payoff(S,optionData,spotData,T1=optionData.T)
	
Where:\n
		S           = Paths of the Underlying.
		optionData  = Datas of the Option.
		spotData  = Datas of the Spot.
		T1 [Optional, default to optionData.T]=Final Time of Spot Simulation (default equals Time to Maturity of the option)

		Payoff      = payoff of the Option.

"""
function payoff(S::AbstractMatrix{num},optionData::AbstractPayoff,spotData::equitySpotData,T1::num2=optionData.T) where{num <: Number,num2 <: Number}
	error("Function used just for documentation")
end
