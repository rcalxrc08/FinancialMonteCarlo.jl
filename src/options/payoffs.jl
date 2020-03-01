####### Payoffs definition

### Spot
include("spot.jl")

### European Payoffs
include("general_european_engine.jl")
include("forward.jl")
include("european_option.jl")
include("binary_european_option.jl")

### Barrier Payoffs
include("barrier_do_option.jl")
include("barrier_di_option.jl")
include("barrier_uo_option.jl")
include("barrier_ui_option.jl")
include("double_barrier_option.jl")

### American Payoffs
include("general_american_engine.jl")
include("american_option.jl")
include("binary_american_option.jl")

### Path Dependent Payoffs
include("general_path_dependent_engine.jl")
include("bermudan_option.jl")
include("general_bermudan_engine.jl")

### Asian Payoffs
include("asian_fixed_strike_option.jl")
include("asian_floating_strike_option.jl")

### Basket Payoffs
include("basket/n_european_option.jl")

### Operations and Strategies
include("operations.jl")
include("operations2.jl")
