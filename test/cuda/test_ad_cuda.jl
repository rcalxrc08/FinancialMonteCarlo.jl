# Automatic differentiation of prices simulated on the GPU (CudaMode), with the backends of the CPU suite: ForwardDiff,
# Enzyme (forward and reverse) and Mooncake (forward and reverse). Run by hand, it needs a GPU and Julia ≤ 1.13.3 for
# Enzyme:
#     julia --project=<env with CUDA, Enzyme, Mooncake, ForwardDiff, Test and this package> test/cuda/test_ad_cuda.jl
# The reference is a central finite difference on the same paths (the seed is reset at every call). The rules of the
# kernels (src/models/adjoints.jl) are for arrays on the host: on the GPU every backend differentiates the generic code,
# and the combinations that do not work yet are marked broken (with the error they throw printed) or skipped.
using Test, FinancialMonteCarlo, CUDA, ForwardDiff, Mooncake, Enzyme
const FMC = FinancialMonteCarlo

const T, K = 1.0, 100.0
const mc = MonteCarloConfiguration(2_000, 10; mode = FMC.CudaMode(), seed = 1)
const x0 = [100.0, 0.2, 0.02, 0.01] # S0, σ, r, d
const models = (BlackScholes = x -> BlackScholesProcess(x[2], Underlying(x[1], x[4])), Bachelier = x -> BachelierProcess(100 * x[2], Underlying(x[1], x[4])), Heston = x -> HestonProcess(0.5, x[2], 0.0, 2.0, -0.7, 0.04, Underlying(x[1], x[4])), Kou = x -> KouProcess(x[2], 3.0, 0.3, 30.0, 20.0, Underlying(x[1], x[4])), LogNormalMixture = x -> LogNormalMixture((x[2], 0.3), (0.5,), Underlying(x[1], x[4])))

# Central differences on the same paths: the payoff is smooth in the parameters along the paths (European call)
function finite_differences(f, x; h = 1e-5)
    return [(f(x .+ h .* (1:length(x) .== i)) - f(x .- h .* (1:length(x) .== i))) / (2h) for i in eachindex(x)]
end

const backends = (ForwardDiff = f -> ForwardDiff.gradient(f, x0), EnzymeForward = f -> collect(Enzyme.gradient(Enzyme.Forward, Enzyme.Const(f), x0)[1]), EnzymeReverse = f -> Enzyme.gradient(Enzyme.Reverse, Enzyme.Const(f), x0)[1], MooncakeReverse = f -> Mooncake.value_and_gradient!!(Mooncake.prepare_gradient_cache(f, x0), f, x0)[2][2], MooncakeForward = function (f)
    cache = Mooncake.prepare_derivative_cache(f, x0)
    return [Mooncake.value_and_derivative!!(cache, (f, Mooncake.zero_tangent(f)), (x0, Float64.(1:4 .== i)))[2] for i = 1:4]
end)

# Backends known not to work on the GPU (CUDA 6.4.1, Enzyme 0.13.219, Mooncake 0.5.70, Julia 1.12.7), for every model
# (#68; ForwardDiff works):
# - EnzymeForward: CUDA's Enzyme extension calls `EnzymeRules.batch_width`, missing from every EnzymeCore, in the
#   batched rules that `Enzyme.gradient(Forward)` reaches (fixed in CUDA 6.4.2, which needs LLVM.jl 10: Enzyme does not
#   allow it yet, EnzymeAD/Enzyme.jl#3741); see also EnzymeAD/Enzyme.jl#4048. One direction at a time (`autodiff` with
#   `Duplicated`) runs, but is silently wrong, primal included: a host-to-device `copyto!` (`CuArray(::Array)`) writes
#   the zero tangent into the primal buffer (`convert(CuPtr, ::Managed)` loses the shadow, fixed in CUDA 6.5,
#   JuliaGPU/CUDA.jl#3332). Do not use Enzyme forward mode on the GPU before CUDA 6.5.
# - EnzymeReverse: no method matching `MixedDuplicated(::Broadcasted{CuArrayStyle...})`: any Float64 scalar in a
#   broadcast over a `Duplicated` CuArray (`y .= s .* x`), EnzymeAD/Enzyme.jl#3451 (fix in progress, #3704), with or
#   without runtime activity.
# - MooncakeReverse: a typeassert on the `IdDict` of a CUDA call (`CoDual{IdDict, NoFData}` instead of `IdDict` fdata).
# - MooncakeForward: building the forward rule crashes the process (illegal instruction in `get_task_tls`), so it is
#   skipped, not run.
const BROKEN = (:EnzymeForward, :EnzymeReverse, :MooncakeReverse)
const SKIPPED = (:MooncakeForward,)

@testset "AD in CudaMode" begin
    @testset "$model" for model in keys(models)
        build = models[model]
        f(x) = pricer(build(x), ZeroRate(x[3]), mc, EuropeanOption(T, K))
        g_ref = finite_differences(f, x0)
        @testset "$backend" for backend in keys(backends)
            if backend in SKIPPED
                @test_skip false
                continue
            end
            g = try
                backends[backend](f)
            catch err
                @info "$model, $backend fails" error = first(sprint(showerror, err), 600)
                nothing
            end
            works = g !== nothing && isapprox(g, g_ref; rtol = 1e-4)
            if backend in BROKEN
                @test_broken works
            else
                @test works
            end
        end
    end
end
