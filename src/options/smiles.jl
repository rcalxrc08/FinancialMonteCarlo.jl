### European Payoffs
include("smiles/european_option_smile.jl")
include("smiles/binary_european_option_smile.jl")

@doc """
    EuropeanOptionSmile(T, K::Vector, isCall::Vector{Bool})

Set of [`EuropeanOption`](@ref)s with the same maturity `T` and strikes `K`, each a call or a put
according to `isCall`. They are priced together on the same paths: [`pricer`](@ref), [`variance`](@ref),
[`confinter`](@ref), [`delta`](@ref) and [`rho`](@ref) return one value per strike. In `CudaMode` the
strikes are computed in one kernel and the results are a `CuArray`, except for the tuples of [`confinter`](@ref).
""" EuropeanOptionSmile

@doc """
    BinaryEuropeanOptionSmile(T, K::Vector, isCall::Vector{Bool})

Set of [`BinaryEuropeanOption`](@ref)s with the same maturity `T` and strikes `K`, each a call or a
put according to `isCall`. They are priced together on the same paths: [`pricer`](@ref), [`variance`](@ref),
[`confinter`](@ref), [`delta`](@ref) and [`rho`](@ref) return one value per strike. In `CudaMode` the
strikes are computed in one kernel and the results are a `CuArray`, except for the tuples of [`confinter`](@ref).
""" BinaryEuropeanOptionSmile
