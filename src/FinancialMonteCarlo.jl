__precompile__()

module FinancialMonteCarlo

using Random # for randn! and related
using LinearAlgebra # for Longstaff Schwartz
using Accessors # setting fields
using ChainRulesCore
include("utils.jl")
include("models.jl")
include("options.jl")
include("io_utils.jl")
include("operations.jl")
include("metrics.jl")
include("output_type.jl")
include("multi_threading.jl")
export pricer, variance, confinter, delta, simulate, payoff, →

end#End Module
