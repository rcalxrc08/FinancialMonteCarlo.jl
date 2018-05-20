## Payoffs
abstract type PayoffMC end
#### Option Data
abstract type OptionData end

abstract type EuropeanPayoff<:PayoffMC end

abstract type AbstractEuropeanOptionData<:OptionData end

####### Payoffs definition

### European Payoffs
include("payoffs/forward.jl")
include("payoffs/european_option.jl")
include("payoffs/binary_european_option.jl")

### Barrier Payoffs
include("payoffs/barrier_option.jl")
include("payoffs/double_barrier_option.jl")

### American Payoffs
include("payoffs/general_american_option.jl")
include("payoffs/american_option.jl")
include("payoffs/binary_american_option.jl")

### Asian Payoffs
include("payoffs/asian_fixed_strike_option.jl")
include("payoffs/asian_floating_strike_option.jl")