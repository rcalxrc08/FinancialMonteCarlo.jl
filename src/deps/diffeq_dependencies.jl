using .DifferentialEquations

include("../models/diff_eq_monte_carlo.jl")

function predict_output_type_zero(::absdiffeqmodel) where {absdiffeqmodel <: MonteCarloDiffEqModel}
    return zero(Float64)
end
