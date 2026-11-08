############### Display Function

import Base.Multimedia.display;
import Base.Multimedia.print;

function print(p::Union{AbstractMonteCarloProcess, AbstractPayoff})
    fldnames = collect(fieldnames(typeof(p)))
    Base.print(typeof(p), "(")
    if (length(fldnames) > 0)
        Base.print(fldnames[1], "=", getfield(p, fldnames[1]))
        popfirst!(fldnames)
        for name in fldnames
            Base.print(",", name, "=", getfield(p, name))
        end
    end
    print(")")
end

function display(p::Union{AbstractMonteCarloProcess, AbstractPayoff})
    print(p)
    println()
end

# Agnostic getters and setters, useful for testing and calibration

"""
    get_parameters(x)

Numeric parameters of a model or of a payoff, flattened into a vector in field order (the
underlying of a model is excluded). The inverse of [`FinancialMonteCarlo.set_parameters!`](@ref).
"""
function get_parameters(x::Union{AbstractPayoff, BaseProcess})
    return parameter_vector(x, parameter_values(x))
end

# The parameter fields of a type: all its fields but the underlying of a model and the tuple of the models of a
# copula. Generated, so that the names are constants and the fields read and set with them are inferred: a loop over
# `fieldnames` gives a Symbol known only at run time, and makes get_parameters and set_parameters! type unstable.
parameter_fields(T::DataType) = Tuple(name for (name, type) in zip(fieldnames(T), fieldtypes(T)) if !(type <: Union{AbstractUnderlying, Tuple}))
@generated parameter_names(::Type{T}) where {T} = QuoteNode(parameter_fields(T))
@generated parameter_values(x) = Expr(:tuple, (:(getfield(x, $(QuoteNode(name)))) for name in parameter_fields(x))...)

# The parameters flattened into one vector of their promoted number type (a Bool `isCall` becomes a number)
parameter_vector(x, ::Tuple{}) = price_type(x)[]
parameter_vector(x, values_::Tuple) = vcat(map(flat_parameter, values_)...)
flat_parameter(x::Number) = x
flat_parameter(x::AbstractArray) = vec(x)

"""
    set_parameters!(model, params)

Return a copy of `model` whose parameters are replaced by the vector `params`, in the order of
[`FinancialMonteCarlo.get_parameters`](@ref). Despite the `!`, `model` is not modified: use the
returned value. Useful to write calibration loops that do not depend on the model.
"""
function set_parameters!(model::XProcess, param::Array{num}) where {XProcess <: BaseProcess, num <: Number}
    names = parameter_names(XProcess)
    N1 = length(names)
    ChainRulesCore.@ignore_derivatives @assert N1 == length(param) "Check number of parameters of the model"
    values_ = ntuple(i -> param[i], Val(N1))
    # One call of the constructor with every field, inferred because `names` is a constant.
    return Accessors.setproperties(model, NamedTuple{names}(values_))
end

function set_parameters!(model::LogNormalMixture, param::Array{num}) where {num <: Number}
    N1 = div(length(param) + 1, 2)
    ChainRulesCore.@ignore_derivatives @assert !(N1 * 2 != length(param) + 1 || N1 < 2) "Check number of parameters of the model"
    eta = param[1:N1]
    lam = param[(N1+1):end]
    return Accessors.setproperties(model, (η = eta, λ = lam))
end

function set_parameters!(model::ShiftedLogNormalMixture, param::Array{num}) where {num <: Number}
    N1 = div(length(param), 2)
    ChainRulesCore.@ignore_derivatives @assert !(N1 * 2 != length(param) || N1 < 2) "Check number of parameters of the model"
    eta = param[1:N1]
    lam = param[(N1+1):(2*N1-1)]
    alfa = param[end]
    return Accessors.setproperties(model, (η = eta, λ = lam, α = alfa))
end
