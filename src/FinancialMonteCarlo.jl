__precompile__()

module FinancialMonteCarlo

	using Requires,Random, LinearAlgebra
	function __init__()
		@require CuArrays = "3a865a2d-5b23-5a0f-bc46-62713ec82fae" include("deps/cuda_dependencies.cujl")
		@require ArrayFire = "b19378d9-d87a-599a-927f-45f220a2c452" include("deps/af_dependencies.cujl")
		@require DualNumbers = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74" include("deps/dual_dependencies.jl")
		@require DifferentialEquations = "0c46a032-eb83-5123-abaf-570d42b7fbaa" include("deps/diffeq_dependencies.jl")
		@require Distributed = "8ba89e20-285c-5b6f-9357-94700520ee1b" include("deps/distributed_dependencies.jl")
	end
	include("utils.jl")
	include("models.jl")
	include("options.jl")
	include("io_utils.jl")
	include("metrics.jl");
	include("output_type.jl")
	include("operations.jl")
	include("multi_threading.jl");
	export
	    pricer,
	    variance,
		confinter,
		delta,
		simulate,
		payoff



end#End Module
