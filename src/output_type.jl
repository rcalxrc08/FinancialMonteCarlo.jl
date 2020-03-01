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

function predict_output_type_zero(mcProcess::FinancialMonteCarlo.BaseProcess)
	model_par=FinancialMonteCarlo.get_parameters(mcProcess);
	zero_out_=sum(model_par)
	
	return zero_out_*zero(zero_out_)*predict_output_type_zero(mcProcess.underlying);

end

function predict_output_type_zero(rfCurve::ZeroRate)
	
	return zero(rfCurve.r);

end

function predict_output_type_zero(und::FinancialMonteCarlo.AbstractUnderlying)
	
	return predict_output_type_zero(und.S0)+predict_output_type_zero(und.d);

end

function predict_output_type_zero(und::FinancialMonteCarlo.Curve)
	
	return zero(eltype(keys_(und)))+zero(eltype(values_(und)));

end

function predict_output_type_zero(rfCurve::FinancialMonteCarlo.AbstractZeroRateCurve)
	
	return predict_output_type_zero(rfCurve.r);

end


function predict_output_type_zero(mcConfig::MonteCarloConfiguration)
	if typeof(mcConfig.parallelMode) <: SerialMode #Move to Cuda section
		return zero(Float64)
	else
		return zero(Float32)
	end
end

function predict_output_type_zero(abstractPayoff::FinancialMonteCarlo.AbstractPayoff)
	opt_par=FinancialMonteCarlo.get_parameters(abstractPayoff);
	zero_out_=sum(opt_par)
	
	return zero_out_*zero(zero_out_);

end

function predict_output_type_zero(abstractPayoff::FinancialMonteCarlo.Spot)
	
	return zero(Int8);

end