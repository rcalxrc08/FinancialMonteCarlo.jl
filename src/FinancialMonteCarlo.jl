__precompile__()

module FinancialMonteCarlo

using Random # for randn! and related
using LinearAlgebra # isposdef and I for the correlation of the copulas
using Accessors # setting fields
import Adapt # moving the arrays of the inputs to the device of the mode
using ChainRulesCore
using MuladdMacro: @muladd # fused multiply-adds in the kernels of the paths
using StaticArrays: SVector, pushfirst # the nodes of the curves: isbits (#79)
using OrderedCollections: OrderedDict # the entries of portfolios and market data sets, in insertion order
include("01_utils/utils.jl")
include("02_models/models.jl")
include("03_options/options.jl")
include("04_operations/io_utils.jl")
include("04_operations/operations.jl")
include("05_metrics/metrics.jl")
include("05_metrics/output_type.jl")
include("06_batches/batches.jl")
include("06_batches/multi_threading.jl")
export pricer, variance, confinter, delta, simulate, payoff, →

end#End Module
