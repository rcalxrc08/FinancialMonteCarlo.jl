import Base: +, -, *, /

# The key of an underlying: its name, or the tuple of the names of the components of a multivariate model or basket
const UnderlyingKey = Union{String, Tuple{Vararg{String}}}

underlying_key(names::Tuple{Vararg{AbstractString}}) = length(names) == 1 ? String(only(names)) : map(String, names)
function underlying_key(name::AbstractString)
    occursin('_', name) || return String(name)
    Base.depwarn("Names joined by `_` (\"$name\") are deprecated: use a tuple of names, $(Tuple(split(name, '_'))).", :→)
    return underlying_key(Tuple(split(name, '_')))
end
underlying_names(key::String) = (key,)
underlying_names(key::Tuple) = key

"""
    MarketDataSet

The models of the underlyings, by name: built with `name → model` (see [`→`](@ref)) and joined with `+`. A multivariate
model is linked to the tuple of the names of its components, e.g. `("ENI", "AAPL") → copula`. The models are kept in
the order they were added; `market[name]` is the model of `name`.
"""
struct MarketDataSet
    models::OrderedDict{UnderlyingKey, AbstractMonteCarloProcess}
end
export MarketDataSet
MarketDataSet() = MarketDataSet(OrderedDict{UnderlyingKey, AbstractMonteCarloProcess}())

"""
    name → model
    name → payoff

Link an underlying to a model, giving a [`MarketDataSet`](@ref), or to a payoff or
[`FinancialMonteCarlo.Position`](@ref), giving a [`FinancialMonteCarlo.Portfolio`](@ref). `name` is a `String` for a single name, and a
tuple of names, one per component, for a multivariate model or a basket payoff: `("ENI", "AAPL")`. Names joined by
`_` (`"ENI_AAPL"`) are deprecated. Type `→` as `\to<TAB>`.
"""
function →(name, model::AbstractMonteCarloProcess)
    key = underlying_key(name)
    key isa String || throw(ArgumentError("a single name model needs a single name, not $key"))
    return MarketDataSet(OrderedDict{UnderlyingKey, AbstractMonteCarloProcess}(key => model))
end
function →(name, model::NDimensionalMonteCarloProcess)
    key = underlying_key(name)
    (key isa Tuple && length(key) == length(model.models)) || throw(ArgumentError("a model of $(length(model.models)) underlyings needs a tuple of $(length(model.models)) names, not $key"))
    return MarketDataSet(OrderedDict{UnderlyingKey, AbstractMonteCarloProcess}(key => model))
end

# Joins market data sets: a name can appear only once, alone or in a multivariate model
function +(x::MarketDataSet, y::MarketDataSet)
    names_x = Set(Iterators.flatten(map(underlying_names, collect(keys(x.models)))))
    for key in keys(y.models)
        any(in(names_x), underlying_names(key)) && throw(ArgumentError("the underlyings of $key are already in the market data set"))
    end
    return MarketDataSet(merge(x.models, y.models))
end

Base.getindex(market::MarketDataSet, name) = market.models[underlying_key(name)]
Base.haskey(market::MarketDataSet, name) = haskey(market.models, underlying_key(name))
Base.iterate(market::MarketDataSet, state...) = iterate(market.models, state...)
Base.length(market::MarketDataSet) = length(market.models)
Base.keys(market::MarketDataSet) = keys(market.models)
Base.values(market::MarketDataSet) = values(market.models)
Base.pairs(market::MarketDataSet) = pairs(market.models)
Base.:(==)(x::MarketDataSet, y::MarketDataSet) = x.models == y.models

number_type(market::MarketDataSet) = number_type(collect(values(market.models)))

function Base.show(io::IO, ::MIME"text/plain", market::MarketDataSet)
    print(io, "MarketDataSet with ", length(market), " underlying", length(market) == 1 ? "" : "s", isempty(market.models) ? "" : ":")
    for (key, model) in market.models
        print(io, "\n  ", key, " → ", nameof(typeof(model)))
    end
end
