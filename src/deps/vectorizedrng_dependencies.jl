using .VectorizedRNG
inner_seed!(tmp_rng::VectorizedRNG.AbstractVRNG,seed)= VectorizedRNG.seed!(seed);