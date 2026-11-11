"""
    Portfolio

The positions on each underlying, by name: built with `name → payoff` or `name → position` (see [`→`](@ref)) and joined
with `+`, which adds the positions on the same name. A basket payoff is linked to the tuple of the names of its
underlyings. Priced against a [`MarketDataSet`](@ref).
"""
struct Portfolio
    positions::OrderedDict{UnderlyingKey, Position}
end
Portfolio() = Portfolio(OrderedDict{UnderlyingKey, Position}())

→(name, position::Position) = Portfolio(OrderedDict{UnderlyingKey, Position}(underlying_key(name) => position))
function →(name, payoff_::SingleNamePayoff)
    key = underlying_key(name)
    key isa String || throw(ArgumentError("a single name payoff needs a single name, not $key"))
    return key → (1 * payoff_)
end
function →(name, payoff_::BasketPayoff)
    key = underlying_key(name)
    key isa Tuple || throw(ArgumentError("a basket payoff needs a tuple of names, one per underlying, not $key"))
    return key → (1 * payoff_)
end

function +(x::Portfolio, y::Portfolio)
    positions = copy(x.positions)
    for (key, position) in y.positions
        positions[key] = haskey(positions, key) ? positions[key] + position : position
    end
    return Portfolio(positions)
end

Base.getindex(portfolio::Portfolio, name) = portfolio.positions[underlying_key(name)]
Base.haskey(portfolio::Portfolio, name) = haskey(portfolio.positions, underlying_key(name))
Base.iterate(portfolio::Portfolio, state...) = iterate(portfolio.positions, state...)
Base.length(portfolio::Portfolio) = length(portfolio.positions)
Base.keys(portfolio::Portfolio) = keys(portfolio.positions)
Base.values(portfolio::Portfolio) = values(portfolio.positions)
Base.pairs(portfolio::Portfolio) = pairs(portfolio.positions)
Base.:(==)(x::Portfolio, y::Portfolio) = x.positions == y.positions

number_type(portfolio::Portfolio) = number_type(collect(values(portfolio.positions)))

function Base.show(io::IO, ::MIME"text/plain", portfolio::Portfolio)
    print(io, "Portfolio with ", length(portfolio), " underlying", length(portfolio) == 1 ? "" : "s", isempty(portfolio.positions) ? "" : ":")
    for (key, position) in portfolio.positions
        print(io, "\n  ", key, " → ", position)
    end
end

# The positions of the portfolio priced on the model of key: the position of a single name, or for a multivariate model
# one position per subset of its components (in the order of compute_indices; undefined where there is none). Nothing if
# the portfolio has no position on these underlyings.
positions_on(portfolio::Portfolio, key::String) = get(portfolio.positions, key, nothing)
function positions_on(portfolio::Portfolio, key::Tuple)
    subsets = compute_indices(length(key))
    out = Vector{Position}(undef, length(subsets))
    found = false
    for (i, subset) in enumerate(subsets)
        name = underlying_key(key[subset])
        if haskey(portfolio.positions, name)
            out[i] = portfolio.positions[name]
            found = true
        end
    end
    return found ? out : nothing
end
