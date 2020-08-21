get_rng_type(x::T) where {T <: AbstractFloat}=zero(T);
get_rng_type(x)=zero(Float64);