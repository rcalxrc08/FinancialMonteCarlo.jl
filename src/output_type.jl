# The number type of the paths and of the prices: the type of the sum of zeros of the number types of the inputs, so that
# dual numbers, Float32 and the other number types of automatic differentiation propagate. It is computed from types only:
# no array is traversed (the strikes of a smile on the GPU included) and nothing is recorded on the tape of ReverseDiff.

"""
    price_type(inputs...)

Number type of the paths and the prices computed from `inputs` (models, curves, payoffs, smiles, methods, numbers,
arrays of them): the type of `zero(T1) + zero(T2) + …` over the number types `Ti` of their parameters, so a dual number
or `Float32` in any input reaches the result. Allocate with it, never with `Float64`.

A new input type needs a method of `FinancialMonteCarlo.number_type` only if it is not a subtype of `BaseProcess{T}`
or `AbstractPayoff{T}`, which carry their number type.
"""
price_type(inputs...) = foldl(promote_number, map(number_type, inputs); init = Int8)

# Int8 is the neutral element: zero(Int8) + zero(T) is a T for every number type of interest
promote_number(::Type{A}, ::Type{B}) where {A <: Number, B <: Number} = typeof(zero(A) + zero(B))

# The number type of one input
number_type(::Type{T}) where {T <: Number} = T
number_type(::T) where {T <: Number} = T
number_type(::AbstractArray{T}) where {T <: Number} = T
# An array of payoffs of one number type, without looking at them
number_type(::AbstractArray{<:AbstractPayoff{num}}) where {num <: Number} = num
number_type(::AbstractArray{Spot}) = Int8
# An array of models or of payoffs of different number types, a Position or a Portfolio: the types of what they hold,
# weights included
number_type(x::Union{AbstractArray, Tuple, AbstractDict}) = foldl(promote_number, (number_type(y) for y in x); init = Int8)
number_type(x::Pair) = price_type(x.first, x.second)

# A method holds no prices: Nsim, Nstep and the seeds do not change their type
number_type(::AbstractMethod) = Int8

number_type(und::AbstractUnderlying) = price_type(und.S0, und.d)
number_type(::CurveType{num1, num2}) where {num1 <: Number, num2 <: Number} = promote_number(num1, num2)
number_type(rfCurve::AbstractZeroRateCurve) = number_type(rfCurve.r)

number_type(::AbstractPayoff{num}) where {num <: Number} = num
# Spot has no parameters: its Float16 must not lower the type of the others
number_type(::Spot) = Int8
number_type(smile::AbstractSmile) = price_type(smile.T, smile.K)

number_type(model::ScalarMonteCarloProcess{num}) where {num <: Number} = price_type(num, model.underlying)
number_type(::Union{AbstractMonteCarloEngine{num}, VectorialMonteCarloProcess{num}}) where {num <: Number} = num
# The paths come from the solver of DifferentialEquations, in Float64 whatever the underlying
number_type(::MonteCarloDiffEqModel) = Float64

##Vectorization "utils"
Base.broadcastable(x::Union{BaseProcess, AbstractZeroRateCurve, AbstractPayoff, AbstractSmile, AbstractMonteCarloConfiguration}) = Ref(x)
