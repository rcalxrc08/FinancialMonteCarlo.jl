
function predict_output_type_zero(x...)::Number
    zero_out_ = sum(y -> predict_output_type_zero(y), x)

    return zero_out_
end

function predict_output_type_zero(x::abstractArray)::Number where {abstractArray <: AbstractArray{T}} where {T}
    zero_out_ = sum(y -> predict_output_type_zero(y), x)

    return zero_out_
end

function predict_output_type_zero(::num)::num where {num <: Number}
    return zero(num)
end

function predict_output_type_zero(und::AbstractUnderlying)::Number
    return predict_output_type_zero(und.S0) + predict_output_type_zero(und.d)
end

function predict_output_type_zero(::CurveType{num1, num2})::Number where {num1 <: Number, num2 <: Number}
    return zero(num1) + zero(num2)
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

function predict_output_type_zero(::baseProcess) where {baseProcess <: AbstractMonteCarloEngine{num}} where {num <: Number}
    return zero(num)
end

function predict_output_type_zero(::baseProcess) where {baseProcess <: VectorialMonteCarloProcess{num}} where {num <: Number}
    return zero(num)
end

##Vectorization "utils"
Base.broadcastable(x::Union{BaseProcess, AbstractZeroRateCurve, AbstractPayoff, AbstractMonteCarloConfiguration}) = Ref(x)
