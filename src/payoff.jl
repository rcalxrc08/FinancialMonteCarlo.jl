## Payoffs
abstract type PayoffMC end
#### Option Data
abstract type OptionData end

### Payoffs definition
include("payoffs/europeanOption.jl")
include("payoffs/asianFloatingStrikeOption.jl")
include("payoffs/asianFixedStrikeOption.jl")
include("payoffs/barrierOption.jl")
include("payoffs/forward.jl")