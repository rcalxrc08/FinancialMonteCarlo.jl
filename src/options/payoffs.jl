####### Payoffs definition

### Spot
include("payoffs/spot.jl")

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
include("payoffs/american_option.jl")
include("payoffs/binary_american_option.jl")

### Path Dependent Payoffs
include("payoffs/bermudan_option.jl")

### Asian Payoffs
include("payoffs/asian_fixed_strike_option.jl")
include("payoffs/asian_floating_strike_option.jl")

### Basket Payoffs
include("payoffs/basket/n_european_option.jl")

### Operations and Strategies
include("operations.jl")
include("operations2.jl")
