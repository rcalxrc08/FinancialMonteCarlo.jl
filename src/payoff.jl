## Payoffs
abstract type PayoffMC end
#### Option Data
abstract type OptionData end

abstract type EuropeanPayoff<:PayoffMC end

abstract type AbstractEuropeanOptionData<:OptionData end

### Payoffs definition
include("payoffs/europeanOption.jl")
include("payoffs/generalAmericanOption.jl")
include("payoffs/americanOption.jl")
include("payoffs/binaryEuropeanOption.jl")
include("payoffs/asianFloatingStrikeOption.jl")
include("payoffs/asianFixedStrikeOption.jl")
include("payoffs/barrierOption.jl")
include("payoffs/forward.jl")