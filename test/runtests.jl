# Test entry point.
#
# Every file runs inside its own module, so files cannot leak variables into each other.
# Groups can be selected by name (case-insensitive), e.g.
#     julia --project=. -e 'using Pkg; Pkg.test(test_args = ["models", "ad"])'
#     julia --project=test test/runtests.jl models
using Test

const TEST_GROUPS = [
    "Structs" => ["structs.jl", "engines.jl"],
    "Models" => ["models/black_scholes.jl", "models/heston.jl", "models/kou.jl", "models/merton.jl", "models/variance_gamma.jl", "models/normal_inverse_gaussian.jl", "models/lognormal_mixture.jl", "models/diffeq.jl"],
    "Methods" => ["methods/sobol.jl", "methods/control_variates.jl", "methods/parallel.jl"],
    "Metrics" => ["metrics/statistics.jl", "metrics/greeks.jl"],
    "AD" => ["ad/dual_numbers.jl", "ad/forward_diff.jl", "ad/reverse_diff.jl", "ad/taylor_series.jl"],
    "Portfolio" => ["portfolio/algebra.jl", "portfolio/multivariate.jl"],
]

const SELECTED = lowercase.(ARGS)

function run_test_file(file)
    mod = @eval Main module $(Symbol(replace(file, r"[/.]" => "_"))) end
    @testset "$file" begin
        Base.include(mod, joinpath(@__DIR__, file))
    end
end

@testset "FinancialMonteCarlo" begin
    for (group, files) in TEST_GROUPS
        isempty(SELECTED) || lowercase(group) in SELECTED || continue
        @testset "$group" begin
            foreach(run_test_file, files)
        end
    end
end
