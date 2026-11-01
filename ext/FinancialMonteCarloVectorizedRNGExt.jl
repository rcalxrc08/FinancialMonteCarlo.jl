module FinancialMonteCarloVectorizedRNGExt

using VectorizedRNG
import FinancialMonteCarlo: batch_rng, inner_seed!

inner_seed!(tmp_rng::VectorizedRNG.AbstractVRNG, seed) = VectorizedRNG.seed!(seed);
batch_rng(::VectorizedRNG.AbstractVRNG) = throw(ArgumentError("MultiThreading needs a generator that can be copied: VectorizedRNG has none, and seeding it resets all of its generators"))

end
