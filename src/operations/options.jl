
# Payoffs are keys of positions: equal when they are the same contract, i.e. the same type (whatever its parameters, so
# that K = 100 and K = 100.0 are one contract) with equal parameters. Structural and allocation free, consistent with
# the two-argument hash that containers use.
Base.hash(x::AbstractPayoff, h::UInt) = hash(parameter_values(x), hash(Base.typename(typeof(x)), h))
Base.isequal(x::AbstractPayoff, y::AbstractPayoff) = Base.typename(typeof(x)) === Base.typename(typeof(y)) && isequal(parameter_values(x), parameter_values(y))
Base.:(==)(x::AbstractPayoff, y::AbstractPayoff) = Base.typename(typeof(x)) === Base.typename(typeof(y)) && parameter_values(x) == parameter_values(y)

# A leg whose weight is exactly zero is dropped. A weight of another number type (a dual number) is kept even when its
# value is zero, since it can carry a derivative: DualNumbers has iszero(dual(0, 1)) == true.
is_null_weight(weight::Union{Integer, AbstractFloat, Rational}) = iszero(weight)
is_null_weight(::Number) = false

"""
    Position

Linear combination of payoffs on the same underlying: payoffs with their weights, in the order they were added.
Positions are built with `+`, `-`, `*` and `/` on payoffs and numbers, e.g. `2 * EuropeanOption(T, K) - Forward(T)`,
or with `Position(payoff => weight, ...)`. The same contract appears once (its weights add up), and a leg whose
weight is zero is removed. A position is priced on one simulation, path by path, so `pricer`, `variance` and
`confinter` of a position are those of its discounted payoff.

A position iterates as `payoff => weight` pairs; `position[payoff]` is the weight of `payoff` (zero if absent).
"""
struct Position{P <: AbstractPayoff, W <: Number}
    # Two vectors, one entry per leg: a position has few legs, so finding a payoff by comparison is cheap, and plain
    # vectors are what every AD back-end (Enzyme included) differentiates through
    payoffs::Vector{P}
    weights::Vector{W}
    function Position(payoffs::Vector{P}, weights::Vector{W}) where {P <: AbstractPayoff, W <: Number}
        ChainRulesCore.@ignore_derivatives length(payoffs) == length(weights) || throw(ArgumentError("one weight per payoff"))
        kept = ChainRulesCore.@ignore_derivatives findall(weight -> !is_null_weight(weight), weights)
        return new{P, W}(payoffs[kept], weights[kept])
    end
end

# A payoff is a position of weight 1
payoff_type(x::AbstractPayoff) = typeof(x)
payoff_type(::Position{P}) where {P} = P
payoff_type(x::Pair) = typeof(x.first)
weight_type(::AbstractPayoff) = Int
weight_type(::Position{P, W}) where {P, W} = W
weight_type(x::Pair) = typeof(x.second)
legs_of(x::AbstractPayoff) = ((x, 1),)
legs_of(x::Position) = zip(x.payoffs, x.weights)
legs_of(x::Pair) = ((x.first, x.second),)

# The sum of payoffs, positions and pairs payoff => weight: a payoff type that holds all of them, the promoted weight
function combine(terms...)
    P = mapreduce(payoff_type, typejoin, terms)
    W = mapreduce(weight_type, promote_type, terms)
    payoffs = P[]
    weights = W[]
    # foreach unrolls the tuple of terms, which can mix payoffs and positions: a for loop over it iterates a Union
    foreach(terms) do term
        for (payoff_, weight) in legs_of(term)
            i = ChainRulesCore.@ignore_derivatives findfirst(isequal(payoff_), payoffs)
            if isnothing(i)
                push!(payoffs, payoff_)
                push!(weights, weight)
            else
                weights[i] += weight
            end
        end
    end
    return Position(payoffs, weights)
end

Position() = Position(AbstractPayoff[], Int[])
Position(legs::Pair{<:AbstractPayoff, <:Number}...) = combine(legs...)

const PositionTerm = Union{AbstractPayoff, Position}
+(x::PositionTerm, y::PositionTerm) = combine(x, y)
*(x::AbstractPayoff, y::Number) = Position([x], [y])
*(x::Position, y::Number) = Position(x.payoffs, x.weights .* y)
*(y::Number, x::PositionTerm) = x * y
/(x::PositionTerm, y::Number) = x * inv(y)
-(x::PositionTerm) = -1 * x
-(x::PositionTerm, y::PositionTerm) = x + (-1 * y)

# A position behaves as the collection of its legs
Base.iterate(position::Position, i::Integer = 1) = i > length(position) ? nothing : (position.payoffs[i] => position.weights[i], i + 1)
Base.length(position::Position) = length(position.payoffs)
Base.isempty(position::Position) = isempty(position.payoffs)
Base.eltype(::Type{Position{P, W}}) where {P, W} = Pair{P, W}
Base.keys(position::Position) = position.payoffs
Base.values(position::Position) = position.weights
Base.pairs(position::Position) = collect(position)
leg_index(position::Position, payoff_::AbstractPayoff) = findfirst(isequal(payoff_), position.payoffs)
Base.haskey(position::Position, payoff_::AbstractPayoff) = !isnothing(leg_index(position, payoff_))
function Base.getindex(position::Position{P, W}, payoff_::AbstractPayoff) where {P, W}
    i = leg_index(position, payoff_)
    return isnothing(i) ? zero(W) : position.weights[i]
end
Base.copy(position::Position) = Position(copy(position.payoffs), copy(position.weights))
# The order of the legs does not matter for equality, nor for the hash
Base.:(==)(x::Position, y::Position) = length(x) == length(y) && all(payoff_ -> haskey(y, payoff_) && x[payoff_] == y[payoff_], x.payoffs)
Base.isequal(x::Position, y::Position) = length(x) == length(y) && all(payoff_ -> haskey(y, payoff_) && isequal(x[payoff_], y[payoff_]), x.payoffs)
Base.hash(x::Position, h::UInt) = hash(reduce(xor, (hash(leg) for leg in x); init = zero(UInt)), hash(:Position, h))
Base.zero(::Union{Position, Type{<:Position}}) = Position()
Base.iszero(position::Position) = isempty(position)

# The largest maturity of the legs, zero for an empty position. The fold starts from a zero of the type of the maturities,
# so that it is type stable (a union of Int8 and Float64 in the fold is an illegal type analysis for Enzyme); when the
# legs are of several types, from an Int8, which does not promote the others.
function maturity(position::Position{P}) where {P}
    T = Base.promote_op(maturity, P)
    return maximum(maturity, position.payoffs; init = isconcretetype(T) ? zero(T) : Int8(0))
end

# The discounted payoff of the position along the paths: the weighted sum of those of its legs
function payoff(S, position::Position, rfCurve::AbstractZeroRateCurve, mcConfig::AbstractMonteCarloConfiguration, T1::Number = maturity(position))
    values_ = position.weights[1] .* payoff(S, position.payoffs[1], rfCurve, mcConfig, T1)
    for i = 2:length(position)
        values_ = values_ .+ position.weights[i] .* payoff(S, position.payoffs[i], rfCurve, mcConfig, T1)
    end
    return values_
end

# The number type of the price: the payoffs and the weights, by their types when the payoffs have one number type
number_type(position::Position{P, W}) where {P, W} = promote_number(legs_number_type(position), W)
legs_number_type(::Position{P}) where {num, P <: AbstractPayoff{num}} = num
legs_number_type(::Position{Spot}) = Int8
legs_number_type(position::Position) = number_type(position.payoffs)

function Base.show(io::IO, position::Position)
    isempty(position) && return print(io, "0")
    for (idx_, (payoff_, weight)) in enumerate(position)
        if weight isa Real
            if idx_ != 1
                print(io, weight > 0 ? " + " : " - ")
            elseif weight < 0
                print(io, "-")
            end
            isone(abs(weight)) || print(io, abs(weight), " * ")
        else
            idx_ != 1 && print(io, " + ")
            print(io, "(", weight, ") * ")
        end
        print(io, payoff_)
    end
end
