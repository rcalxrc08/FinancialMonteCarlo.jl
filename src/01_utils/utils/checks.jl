# Validation of the inputs (#71): explicit exceptions, never @assert, which throws a generic AssertionError (not to be told
# apart from a bug) and may be disabled. DomainError for a value out of its domain, DimensionMismatch for lengths and sizes,
# ArgumentError for the rest. The throw is out of line, as in FinancialFFT and PatternVectors (Enzyme aborted on an inlined
# constructor error, EnzymeAD/Enzyme.jl#4053).
@noinline throw_domain_error(value, msg) = throw(DomainError(value, msg))
@noinline throw_dimension_mismatch(msg) = throw(DimensionMismatch(msg))
@noinline throw_argument_error(msg) = throw(ArgumentError(msg))
check_domain(cond::Bool, value, msg::AbstractString) = cond ? nothing : throw_domain_error(value, msg)
check_dimension(cond::Bool, msg::AbstractString) = cond ? nothing : throw_dimension_mismatch(msg)
check_argument(cond::Bool, msg::AbstractString) = cond ? nothing : throw_argument_error(msg)

# The validations are named functions marked @non_differentiable, called directly with the raw inputs, one per
# constructor (#85): a ChainRulesCore.@ignore_derivatives block builds a closure whose captures Zygote traces, ~11 μs each
# with Zygote 0.7.18 (FluxML/Zygote.jl#1753), ~2.5 μs with 0.7.19, against nothing measurable for a @non_differentiable
# function; FinancialCalibrator rebuilds the model at every iteration. Enzyme and Mooncake go through them as through any
# function with no active output. The conditions are computed inside the functions, so the comparisons are not traced
# either. check_domain, check_dimension and check_argument are non-differentiable too, for the callers of one check.
check_maturity(T) = check_domain(T > 0, T, "Time to Maturity must be positive")
check_volatility(σ) = check_domain(σ > 0, σ, "Volatility must be positive")
check_spot(S0) = check_domain(S0 >= zero(S0), S0, "Underlying starting value must be positive")

# A zero of the promotion of number types, the numtype of a parametric struct: the arithmetic of the inputs (zero(A) +
# zero(B)), not promote_type, so that every number type (nested duals, tracked reals) gets the type its operations give.
zero_of_types(types::Type...) = +(map(zero, types)...)

ChainRulesCore.@non_differentiable check_domain(::Any, ::Any, ::Any)
ChainRulesCore.@non_differentiable check_dimension(::Any, ::Any)
ChainRulesCore.@non_differentiable check_argument(::Any, ::Any)
ChainRulesCore.@non_differentiable check_maturity(::Any)
ChainRulesCore.@non_differentiable check_volatility(::Any)
ChainRulesCore.@non_differentiable check_spot(::Any)
ChainRulesCore.@non_differentiable zero_of_types(::Any...)
