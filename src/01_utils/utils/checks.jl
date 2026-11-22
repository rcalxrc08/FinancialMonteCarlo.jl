# Validation of the inputs (#71): explicit exceptions, never @assert, which throws a generic AssertionError (not to be told
# apart from a bug) and may be disabled. DomainError for a value out of its domain, DimensionMismatch for lengths and sizes,
# ArgumentError for the rest. The callers wrap them in ChainRulesCore.@ignore_derivatives, so AD never sees the branch.
check_domain(cond::Bool, value, msg::AbstractString) = cond ? nothing : throw(DomainError(value, msg))
check_dimension(cond::Bool, msg::AbstractString) = cond ? nothing : throw(DimensionMismatch(msg))
check_argument(cond::Bool, msg::AbstractString) = cond ? nothing : throw(ArgumentError(msg))
