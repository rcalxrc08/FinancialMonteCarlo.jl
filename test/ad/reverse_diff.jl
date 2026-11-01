include("../testutils.jl")
using ForwardDiff, ReverseDiff

@testset "ReverseDiff" begin
    mc_tiny = MonteCarloConfiguration(100, Nstep_small)
    x0 = [sigma, S0, r, d, T]
    builders = (
        fwd = τ -> Forward(τ),
        eu = τ -> EuropeanOption(τ, K),
        am = τ -> AmericanOption(τ, K),
        barrier = τ -> BarrierOptionDownOut(τ, K, D),
        asian_float = τ -> AsianFloatingStrikeOption(τ),
        asian_fixed = τ -> AsianFixedStrikeOption(τ, K),
    )
    @testset "gradient: $name" for (name, build) in pairs(builders)
        f(x) = pricer(BlackScholesProcess(x[1], Underlying(x[2], x[4])), ZeroRate(x[3]), mc_tiny, build(x[5]))
        g = ReverseDiff.gradient(f, x0)
        @test length(g) == length(x0)
        @test all(isfinite, g)
    end

    @testset "gradient: smile" begin
        smile(τ, k) = FinancialMonteCarlo.EuropeanOptionSmile(τ, [k, 1.1 * k], [true, false])
        f(x) = sum(pricer(BlackScholesProcess(x[1], Underlying(x[2], x[4])), ZeroRate(x[3]), mc_tiny, smile(x[5], K)))
        @test ReverseDiff.gradient(f, x0) ≈ ForwardDiff.gradient(f, x0)
    end
end
