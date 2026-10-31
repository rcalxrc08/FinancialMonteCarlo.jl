module FinancialMonteCarloMooncakeExt

using FinancialMonteCarlo, Mooncake
using Mooncake: @zero_derivative, MinimalCtx

# Reseeding the rng ends in a foreigncall for the MersenneTwister (the default rng), and the seed carries no derivative.
@zero_derivative MinimalCtx Tuple{typeof(FinancialMonteCarlo.set_seed!), Any}

end
