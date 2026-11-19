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
parameter_fields(T::DataType) = Tuple(name for (name, type) in zip(fieldnames(T), fieldtypes(T)) if !(type <: Union{AbstractUnderlying, Tuple{Vararg{BaseProcess}}}))
@generated parameter_names(::Type{T}) where {T} = QuoteNode(parameter_fields(T))
@generated parameter_values(x) = Expr(:tuple, (:(getfield(x, $(QuoteNode(name)))) for name in parameter_fields(x))...)

# The parameters flattened into one vector of their promoted number type (a Bool `isCall` becomes a number)
parameter_vector(x, ::Tuple{}) = price_type(x)[]
parameter_vector(x, values_::Tuple) = vcat(map(flat_parameter, values_)...)
flat_parameter(x::Number) = x
flat_parameter(x::AbstractArray) = vec(x)
flat_parameter(x::Tuple) = collect(x)

"""
    set_parameters!(model, params)

Return a copy of `model` whose parameters are replaced by the vector `params`, in the order of
[`FinancialMonteCarlo.get_parameters`](@ref). Despite the `!`, `model` is not modified: use the
returned value. Useful to write calibration loops that do not depend on the model.
"""
function set_parameters!(model::XProcess, param::AbstractVector{num}) where {XProcess <: BaseProcess, num <: Number}
    names = parameter_names(XProcess)
    N1 = length(names)
    ChainRulesCore.@ignore_derivatives check_dimension(N1 == length(param), "Check number of parameters of the model")
    values_ = ntuple(i -> param[i], Val(N1))
    # One call of the constructor with every field, inferred because `names` is a constant.
    return Accessors.setproperties(model, NamedTuple{names}(values_))
end

function set_parameters!(model::LogNormalMixture{N, M}, param::AbstractVector{num}) where {N, M, num <: Number}
    ChainRulesCore.@ignore_derivatives check_dimension(length(param) == N + M, "Check number of parameters of the model")
    return LogNormalMixture(ntuple(i -> param[i], Val(N)), ntuple(i -> param[N+i], Val(M)), model.underlying)
end

# The parameters of the model, then the shift (a method of its own: through the generic one, get_parameters would call
# itself on the inner model, and inference gives up on the recursion)
get_parameters(model::ShiftedProcess) = vcat(get_parameters(model.model), model.α)
function set_parameters!(model::ShiftedProcess, param::AbstractVector{num}) where {num <: Number}
    ChainRulesCore.@ignore_derivatives check_dimension(length(param) == length(get_parameters(model)), "Check number of parameters of the model")
    return ShiftedProcess(set_parameters!(model.model, param[1:(end-1)]), param[end], model.underlying)
end

# Reverse rule of set_parameters! (#80): Zygote through the constructors and the reads of `param` costs 0.1-0.2 ms. The
# pullback reads the tangent of the returned model in the order of get_parameters, and gives the tangent of the fields
# that are not parameters (the underlying, which set_parameters! copies) to the input model.
tangent_number(x) = x
tangent_number(::AbstractZero) = false
tangent_numbers(t, N::Integer) = [tangent_number(t[i]) for i = 1:N]
tangent_numbers(::AbstractZero, N::Integer) = fill(false, N)
tangent_field(Δ, name::Symbol) = getproperty(Δ, name)
tangent_field(Δ::AbstractZero, ::Symbol) = Δ
# generated as parameter_values, so that the fields of the tangent are read with constant names (inferred)
@generated parameter_tangent_values(Δ, ::Type{T}) where {T} = Expr(:tuple, (:(tangent_number(tangent_field(Δ, $(QuoteNode(name))))) for name in parameter_fields(T))...)
parameter_tangent(model::BaseProcess, Δ) = [parameter_tangent_values(Δ, typeof(model))...]
parameter_tangent(::LogNormalMixture{N, M}, Δ) where {N, M} = vcat(tangent_numbers(tangent_field(Δ, :η), N), tangent_numbers(tangent_field(Δ, :λ), M))
parameter_tangent(model::ShiftedProcess, Δ) = vcat(parameter_tangent(model.model, tangent_field(Δ, :model)), tangent_number(tangent_field(Δ, :α)))
underlying_tangent(::BaseProcess, Δ::AbstractZero) = Δ
underlying_tangent(model::BaseProcess, Δ) = Tangent{typeof(model)}(; underlying = tangent_field(Δ, :underlying))
underlying_tangent(::ShiftedProcess, Δ::AbstractZero) = Δ
underlying_tangent(model::ShiftedProcess, Δ) = Tangent{typeof(model)}(; model = underlying_tangent(model.model, tangent_field(Δ, :model)), underlying = tangent_field(Δ, :underlying))
function ChainRulesCore.rrule(::typeof(set_parameters!), model::BaseProcess, param::AbstractVector{<:Number})
    project = ProjectTo(param)
    function set_parameters_pullback(Δ)
        Δ_ = unthunk(Δ)
        return NoTangent(), underlying_tangent(model, Δ_), project(parameter_tangent(model, Δ_))
    end
    return set_parameters!(model, param), set_parameters_pullback
end
