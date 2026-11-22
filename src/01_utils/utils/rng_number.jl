get_rng_type(::T) where {T <: AbstractFloat} = zero(T);
get_rng_type(::Any) = zero(Float64);
