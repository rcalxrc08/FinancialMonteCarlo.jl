function predict_output_type(x...)
	
	zero_out_=predict_output_type_zero(x...)

	return typeof(zero_out_);
end

function predict_output_type_zero(x...)
	
	zero_out_=sum(y->predict_output_type_zero(y),x)

	return zero_out_;
end

function predict_output_type_zero(x::abstractArray) where {abstractArray <: AbstractArray{T}} where T
	
	zero_out_=sum(y->predict_output_type_zero(y),x)

	return zero_out_;
end

function predict_output_type_zero(x::num)::num where {num <: Number}
	return zero(x);
end

function predict_output_type_zero(mcProcess::BaseProcess)
	model_par=get_parameters(mcProcess);
	zero_out_=sum(model_par)
	
	return zero_out_*zero(zero_out_)*predict_output_type_zero(mcProcess.underlying);

end

function predict_output_type_zero(mcProcess::VectorialMonteCarloProcess)
	zero_out_=sum(y->predict_output_type_zero(y),mcProcess.models)
	model_par=sum(get_parameters(mcProcess));
	return zero(zero_out_)*zero(model_par)
	
end

function predict_output_type_zero(und::AbstractUnderlying)
	
	return predict_output_type_zero(und.S0)+predict_output_type_zero(und.d);

end

function predict_output_type_zero(und::Curve)
	
	return zero(eltype(keys_(und)))+zero(eltype(values_(und)));

end

function predict_output_type_zero(rfCurve::AbstractZeroRateCurve)
	
	return predict_output_type_zero(rfCurve.r);

end


function predict_output_type_zero(mcConfig::MonteCarloConfiguration)
	if typeof(mcConfig.parallelMode) <: SerialMode #Move to Cuda section
		return zero(Float64)
	else
		return zero(Float32)
	end
end

function predict_output_type_zero(abstractPayoff::AbstractPayoff)
	opt_par=get_parameters(abstractPayoff);
	zero_out_=sum(opt_par)
	
	return zero_out_*zero(zero_out_);

end

function predict_output_type_zero(abstractPayoff::Spot)
	
	return zero(Int8);

end