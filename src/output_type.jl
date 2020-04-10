function predict_output_type(x...)::Type
	
	zero_out_=predict_output_type_zero(x...)

	return typeof(zero_out_);
end

function predict_output_type_zero_(x...)::Number
	
	zero_out_=sum(y->predict_output_type_zero(y),x)

	return zero_out_;
end

function predict_output_type_zero(x::abstractArray)::Number where {abstractArray <: AbstractArray{T}} where T
	
	zero_out_=sum(y->predict_output_type_zero(y),x)

	return zero_out_;
end

function predict_output_type_zero(x::num)::num where {num <: Number}
	return zero(x);
end

#function predict_output_type_zero(mcProcess::BaseProcess)::Number
#	model_par=get_parameters(mcProcess);
#	zero_out_=sum(model_par)
#	
#	return zero_out_*zero(zero_out_)*predict_output_type_zero(mcProcess.underlying);
#
#end

# function predict_output_type_zero(mcProcess::AbstractMonteCarloEngine)::Number
	# return sum(zero.(fieldtypes(typeof(mcProcess))))
# end

# function predict_output_type_zero(mcProcess::VectorialMonteCarloProcess)::Number
	# zero_out_=sum(y->predict_output_type_zero(y),mcProcess.models)
	# model_par=sum(get_parameters(mcProcess));
	# return zero(zero_out_)*zero(model_par)
	
# end

function predict_output_type_zero(und::AbstractUnderlying)::Number
	
	return predict_output_type_zero(und.S0)+predict_output_type_zero(und.d);

end

function predict_output_type_zero(und::Curve)::Number
	
	return zero(eltype(keys_(und)))+zero(eltype(values_(und)));

end

function predict_output_type_zero(rfCurve::AbstractZeroRateCurve)::Number
	
	return predict_output_type_zero(rfCurve.r);

end


function predict_output_type_zero(mcConfig::MonteCarloConfiguration)::Int8
	zero(Int8)
end

#function predict_output_type_zero(abstractPayoff::AbstractPayoff)::Number
#	opt_par=get_parameters(abstractPayoff);
#	zero_out_=sum(opt_par)
#	
#	return zero_out_*zero(zero_out_);
#
#end

function predict_output_type_zero(abstractPayoff::Spot)::Int8
	
	return zero(Int8);

end


function predict_output_type_zero(x::basePayoff) where {basePayoff <: AbstractPayoff{num}} where {num <: Number}
	return zero(num);
end

function predict_output_type_zero(x::baseProcess) where {baseProcess <: ScalarMonteCarloProcess{num}} where {num <: Number}
	return zero(num)+predict_output_type_zero(x.underlying);
end

function predict_output_type_zero(x::baseProcess) where {baseProcess <: VectorialMonteCarloProcess{num}} where {num <: Number}
	return zero(num);
end

##Vectorization "utils"
import Base.iterate
iterate(x::Union{BaseProcess,AbstractZeroRateCurve,AbstractPayoff,AbstractMonteCarloConfiguration}) where {proc <: BaseProcess}=(x,nothing)
iterate(x::Union{BaseProcess,AbstractZeroRateCurve,AbstractPayoff,AbstractMonteCarloConfiguration},p::Nothing)=nothing

import Base.length
length(x::Union{BaseProcess,AbstractZeroRateCurve,AbstractPayoff,AbstractMonteCarloConfiguration})=1