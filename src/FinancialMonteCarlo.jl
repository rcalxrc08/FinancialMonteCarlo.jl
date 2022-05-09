__precompile__()

module FinancialMonteCarlo

using Requires # for conditional dependencies
using Random # for randn! and related
using LinearAlgebra # for Longstaff Schwartz
function __init__()
    @require CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba" include("deps/cuda_dependencies.cujl")
    @require DifferentialEquations = "0c46a032-eb83-5123-abaf-570d42b7fbaa" include("deps/diffeq_dependencies.jl")
    @require VectorizedRNG = "33b4df10-0173-11e9-2a0c-851a7edac40e" include("deps/vectorizedrng_dependencies.jl")
    @require Distributed = "8ba89e20-285c-5b6f-9357-94700520ee1b" include("deps/distributed_dependencies.jl")
    @require TaylorSeries = "6aa5eb33-94cf-58f4-a9d0-e4b2c4fc25ea" include("deps/taylorseries_dependencies.jl")
end
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
