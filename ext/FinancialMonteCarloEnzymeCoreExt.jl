module FinancialMonteCarloEnzymeCoreExt

using FinancialMonteCarlo
using EnzymeCore: EnzymeRules

# The seed and the random draws carry no derivative. Reseeding a generator can end in a foreigncall, and Enzyme cannot
# analyse the scalar draws of the MersenneTwister (its state is read as floats and as integers), so they are inactive:
# they run in the forward pass and their results are constants.
EnzymeRules.inactive(::typeof(FinancialMonteCarlo.set_seed!), args...) = nothing
EnzymeRules.inactive(::typeof(FinancialMonteCarlo.draw_uniform), args...) = nothing
EnzymeRules.inactive(::typeof(FinancialMonteCarlo.draw_exponential), args...) = nothing
EnzymeRules.inactive(::typeof(FinancialMonteCarlo.draw_normal), args...) = nothing

end
