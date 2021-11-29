function predict_output_type(x...)::Type
    zero_out_ = predict_output_type_zero(x...)

    return typeof(zero_out_)
end

function predict_output_type_zero_(x...)::Number
    zero_out_ = sum(y -> predict_output_type_zero(y), x)

    return zero_out_
end

function predict_output_type_zero(x::abstractArray)::Number where {abstractArray <: AbstractArray{T}} where {T}
    zero_out_ = sum(y -> predict_output_type_zero(y), x)

    return zero_out_
end

function predict_output_type_zero(x::num)::num where {num <: Number}
    return zero(x)
end

function predict_output_type_zero(und::AbstractUnderlying)::Number
    return predict_output_type_zero(und.S0) + predict_output_type_zero(und.d)
end

function predict_output_type_zero(und::CurveType)::Number
    return zero(eltype(keys_(und))) + zero(eltype(values_(und)))
end

function predict_output_type_zero(rfCurve::AbstractZeroRateCurve)::Number
    return predict_output_type_zero(rfCurve.r)
end

function predict_output_type_zero(::Any)::Int8
    zero(Int8)
end

function predict_output_type_zero(::Spot)::Int8
    return zero(Int8)
end

function predict_output_type_zero(::basePayoff) where {basePayoff <: AbstractPayoff{num}} where {num <: Number}
    return zero(num)
end

function predict_output_type_zero(x::baseProcess) where {baseProcess <: ScalarMonteCarloProcess{num}} where {num <: Number}
    return zero(num) + predict_output_type_zero(x.underlying)
end

function predict_output_type_zero(::baseProcess) where {baseProcess <: VectorialMonteCarloProcess{num}} where {num <: Number}
    return zero(num)
end

##Vectorization "utils"
import Base.iterate
iterate(x::Union{BaseProcess, AbstractZeroRateCurve, AbstractPayoff, AbstractMonteCarloConfiguration}) = (x, nothing)
iterate(::Union{BaseProcess, AbstractZeroRateCurve, AbstractPayoff, AbstractMonteCarloConfiguration}, p::Nothing) = nothing

import Base.length
length(x::Union{BaseProcess, AbstractZeroRateCurve, AbstractPayoff, AbstractMonteCarloConfiguration}) = 1
