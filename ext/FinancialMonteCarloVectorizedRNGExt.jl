module FinancialMonteCarloVectorizedRNGExt

using VectorizedRNG
import FinancialMonteCarlo: inner_seed!

inner_seed!(tmp_rng::VectorizedRNG.AbstractVRNG, seed) = VectorizedRNG.seed!(seed);

end
