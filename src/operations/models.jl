import Base.+;

"""
    MarketDataSet

Dictionary from underlying names to models, built with `name → model` and joined with `+`.
The name of a multivariate model joins the names of its components with `_`, e.g. `"ENI_AAPL"`.
"""
const MarketDataSet = Dict{String, AbstractMonteCarloProcess};
export MarketDataSet;

"""
    name → model
    name → payoff

Link the underlying `name` to a model, giving a [`MarketDataSet`](@ref), or to a payoff or
[`FinancialMonteCarlo.Position`](@ref), giving a [`FinancialMonteCarlo.Portfolio`](@ref).
Single name keys cannot contain `_`; multivariate models and basket payoffs need one name per
component, joined by `_`. Type `→` as `\\to<TAB>`.
"""
function →(x::String, y::AbstractMonteCarloProcess)
    out = MarketDataSet(x => y)
    return out
end

#Link an multiname underlying to a multi variate model
function →(x::String, y::NDimensionalMonteCarloProcess)
    sep = findfirst("_", x)
    ChainRulesCore.@ignore_derivatives @assert !isnothing(sep) "Nvariate process must follow the format: INDEX1_INDEX2"
    len_ = length(y.models)
    sep = split(x, "_")
    ChainRulesCore.@ignore_derivatives @assert len_ == length(sep)
    out = MarketDataSet(x => y)

    return out
end

#Joins two different MarketDataSet s
function +(x::MarketDataSet, y::MarketDataSet)
    out = copy(x)
    y_keys = keys(y)
    for y_key in y_keys
        ChainRulesCore.@ignore_derivatives @assert !(haskey(out, y_key) || any(out_el -> any(out_el_sub -> out_el_sub == y_key, split(out_el, "_")), keys(out)) || any(y_key_el -> haskey(out, y_key_el), split(y_key, "_"))) "check keys"
        out[y_key] = y[y_key]
    end
    return out
end
