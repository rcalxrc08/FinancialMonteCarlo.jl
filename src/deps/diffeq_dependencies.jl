using .DifferentialEquations

include("../models/diff_eq_monte_carlo.jl")

function predict_output_type_zero(x::absdiffeqmodel) where {absdiffeqmodel <: MonteCarloDiffEqModel}
	
	return zero(Float64);
end