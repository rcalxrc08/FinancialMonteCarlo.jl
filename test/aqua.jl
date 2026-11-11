include("testutils.jl")
using Aqua

# Method ambiguities, type piracy, unbound type parameters, undefined exports, stale dependencies and compat bounds
@testset "Aqua" begin
    # persistent_tasks precompiles the package in a fresh environment, minutes for a check that only matters for packages
    # starting tasks at load
    Aqua.test_all(FinancialMonteCarlo; persistent_tasks = false)
end
