### European Payoffs
include("smiles/european_option_smile.jl")
include("smiles/binary_european_option_smile.jl")

@doc """
    EuropeanOptionSmile(T, K::Vector, isCall::Vector{Bool})

Set of [`EuropeanOption`](@ref)s with the same maturity `T` and strikes `K`, each a call or a put
according to `isCall`. [`pricer`](@ref) prices them all on the same paths and returns the vector of prices.
""" EuropeanOptionSmile

@doc """
    BinaryEuropeanOptionSmile(T, K::Vector, isCall::Vector{Bool})

Set of [`BinaryEuropeanOption`](@ref)s with the same maturity `T` and strikes `K`, each a call or a
put according to `isCall`. [`pricer`](@ref) prices them all on the same paths and returns the vector of prices.
""" BinaryEuropeanOptionSmile
