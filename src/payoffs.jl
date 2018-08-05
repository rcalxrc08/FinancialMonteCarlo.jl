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