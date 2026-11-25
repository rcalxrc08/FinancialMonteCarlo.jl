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

The market of a portfolio: the underlyings by name (`name → underlying`) and the models by key (`key → model`), see
[`→`](@ref), joined with `+`. A single name model is linked to its name, a multivariate model to the tuple of the names
of its components, e.g. `("ENI", "AAPL") → copula`, whose underlyings are those of the names. `name → (model,
underlying)` links both at once. Models and underlyings are kept in the order they were added; `market[key]` is the
model of `key`, `underlying(market, key)` its underlying (a tuple for a multivariate model).
"""
struct MarketDataSet
    underlyings::OrderedDict{String, AbstractUnderlying}
    models::OrderedDict{UnderlyingKey, AbstractMonteCarloProcess}
end
export MarketDataSet
MarketDataSet() = MarketDataSet(OrderedDict{String, AbstractUnderlying}(), OrderedDict{UnderlyingKey, AbstractMonteCarloProcess}())
model_set(key, model) = MarketDataSet(OrderedDict{String, AbstractUnderlying}(), OrderedDict{UnderlyingKey, AbstractMonteCarloProcess}(key => model))

"""
    name → underlying
    name → model
    name → (model, underlying)
    name → payoff

Link a name to its [`Underlying`](@ref) or to a model, giving a [`MarketDataSet`](@ref), or to a payoff or
[`FinancialMonteCarlo.Position`](@ref), giving a [`FinancialMonteCarlo.Portfolio`](@ref). `name` is a `String` for a single name, and a
tuple of names, one per component, for a multivariate model or a basket payoff: `("ENI", "AAPL")`. Names joined by
`_` (`"ENI_AAPL"`) are deprecated. Type `→` as `\to<TAB>`.
"""
function →(name, model::AbstractMonteCarloProcess)
    key = underlying_key(name)
    key isa String || throw(ArgumentError("a single name model needs a single name, not $key"))
    return model_set(key, model)
end
function →(name, und::AbstractUnderlying)
    key = underlying_key(name)
    key isa String || throw(ArgumentError("an underlying needs a single name, not $key"))
    return MarketDataSet(OrderedDict{String, AbstractUnderlying}(key => und), OrderedDict{UnderlyingKey, AbstractMonteCarloProcess}())
end
# A model and its underlying (a tuple of them, one per name, for a multivariate model) at once
function →(name, (model, und)::Tuple{AbstractMonteCarloProcess, Union{AbstractUnderlying, UnderlyingTuple}})
    names = underlying_names(underlying_key(name))
    unds = und isa AbstractUnderlying ? (und,) : und
    length(unds) == length(names) || throw(ArgumentError("one underlying per name: $(length(names)) names, $(length(unds)) underlyings"))
    return (name → model) + reduce(+, (n → u for (n, u) in zip(names, unds)))
end
function →(name, model::NDimensionalMonteCarloProcess)
    key = underlying_key(name)
    (key isa Tuple && length(key) == length(model.models)) || throw(ArgumentError("a model of $(length(model.models)) underlyings needs a tuple of $(length(model.models)) names, not $key"))
    return model_set(key, model)
end

# Joins market data sets: a name can have one underlying, and appear in one model only, alone or in a multivariate one
function +(x::MarketDataSet, y::MarketDataSet)
    names_x = Set(Iterators.flatten(map(underlying_names, collect(keys(x.models)))))
    for key in keys(y.models)
        any(in(names_x), underlying_names(key)) && throw(ArgumentError("the underlyings of $key already have a model in the market data set"))
    end
    for name in keys(y.underlyings)
        haskey(x.underlyings, name) && throw(ArgumentError("$name already has an underlying in the market data set"))
    end
    return MarketDataSet(merge(x.underlyings, y.underlyings), merge(x.models, y.models))
end

"""
    underlying(market, key)

The [`Underlying`](@ref) of the name `key` in the [`MarketDataSet`](@ref) `market`, or the tuple of the underlyings of
the names of a multivariate key, in the order of the names.
"""
underlying(market::MarketDataSet, key) = market_underlying(market, underlying_key(key))
function market_underlying(market::MarketDataSet, key::String)
    haskey(market.underlyings, key) || throw(ArgumentError("no underlying for $key in the market data set: add it with $key → Underlying(S0, d)"))
    return market.underlyings[key]
end
market_underlying(market::MarketDataSet, key::Tuple) = map(name -> market_underlying(market, name), key)

Base.getindex(market::MarketDataSet, name) = market.models[underlying_key(name)]
Base.haskey(market::MarketDataSet, name) = haskey(market.models, underlying_key(name))
Base.iterate(market::MarketDataSet, state...) = iterate(market.models, state...)
Base.length(market::MarketDataSet) = length(market.models)
Base.keys(market::MarketDataSet) = keys(market.models)
Base.values(market::MarketDataSet) = values(market.models)
Base.pairs(market::MarketDataSet) = pairs(market.models)
Base.:(==)(x::MarketDataSet, y::MarketDataSet) = x.models == y.models

number_type(market::MarketDataSet) = price_type(collect(values(market.models)), collect(values(market.underlyings)))

function Base.show(io::IO, ::MIME"text/plain", market::MarketDataSet)
    print(io, "MarketDataSet with ", length(market), " underlying", length(market) == 1 ? "" : "s", isempty(market.models) ? "" : ":")
    for (key, model) in market.models
        print(io, "\n  ", key, " → ", nameof(typeof(model)))
    end
    for (name, und) in market.underlyings
        print(io, "\n  ", name, " → ", und)
    end
end
