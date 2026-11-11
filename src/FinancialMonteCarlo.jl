__precompile__()

module FinancialMonteCarlo

using Random # for randn! and related
using LinearAlgebra # isposdef and I for the correlation of the copulas
using Accessors # setting fields
import Adapt # moving the arrays of the inputs to the device of the mode
using ChainRulesCore
using OrderedCollections: OrderedDict # the entries of portfolios and market data sets, in insertion order
include("utils.jl")
include("models.jl")
include("options.jl")
include("io_utils.jl")
include("operations.jl")
include("metrics.jl")
include("output_type.jl")
include("batches.jl")
include("multi_threading.jl")
include("deprecated.jl")
export pricer, variance, confinter, delta, simulate, payoff, →

end#End Module
