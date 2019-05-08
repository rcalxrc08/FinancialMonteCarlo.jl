using .DifferentialEquations

include("models/diff_eq_monte_carlo.jl")

pricer_macro(MonteCarloProblem)
pricer_macro_array(MonteCarloProblem)

variance_macro(MonteCarloProblem)
variance_macro_array(MonteCarloProblem)

confinter_macro(MonteCarloProblem)
confinter_macro_array(MonteCarloProblem)