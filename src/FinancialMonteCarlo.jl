__precompile__()

module FinancialMonteCarlo

	using Requires,Random, LinearAlgebra
	function __init__()
		@require CuArrays = "3a865a2d-5b23-5a0f-bc46-62713ec82fae" include("cuda_dependencies.cujl")
		@require ArrayFire = "b19378d9-d87a-599a-927f-45f220a2c452" include("af_dependencies.cujl")
		@require DualNumbers = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74" include("dual_dependencies.jl")
		@require DifferentialEquations = "0c46a032-eb83-5123-abaf-570d42b7fbaa" include("diffeq_dependencies.jl")
	end
	include("utils.jl")
	include("payoffs.jl")
	include("models.jl")
	include("metrics.jl");
	include("multi_threading.jl");
	include("multi_process.jl");
	export
	    pricer,
	    variance,
		confinter,
		delta,
		simulate,
		payoff



end#End Module
