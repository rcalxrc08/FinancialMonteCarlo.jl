include("../testutils.jl")
using Mooncake, ForwardDiff, Random

@testset "Mooncake" begin
    mc_xo = MonteCarloConfiguration(100, Nstep_small, FinancialMonteCarlo.StandardMC(), FinancialMonteCarlo.SerialMode(0, Xoshiro()))

    @testset "Kou European call: gradient matches ForwardDiff" begin
        # x = [S0, σ, λ, p, λ₊, λ₋, r, d, T]
        x0 = [S0, sigma, 5.0, 0.3, 30.0, 20.0, r, d, T]
        f(x) = pricer(KouProcess(x[2], x[3], x[4], x[5], x[6], Underlying(x[1], x[8])), ZeroRate(x[7]), mc_xo, EuropeanOption(x[9], K))
        cache = Mooncake.prepare_gradient_cache(f, x0)
        price, (_, g) = Mooncake.value_and_gradient!!(cache, f, x0)
        @test price ≈ f(x0)
        @test g ≈ ForwardDiff.gradient(f, x0) rtol = 1e-8
        @test 0 < g[1] < 1 # delta
        @test g[2] > 0 # vega
    end

    @testset "default rng: seeding has a rule, refills go through" begin
        # 2000 paths draw past the 1002 numbers cached by the MersenneTwister.
        x0 = [S0, sigma, r, d]
        mc_default = MonteCarloConfiguration(2000, Nstep_small)
        f(x) = pricer(BlackScholesProcess(x[2], Underlying(x[1], x[4])), ZeroRate(x[3]), mc_default, EuropeanOption(T, K))
        cache = Mooncake.prepare_gradient_cache(f, x0)
        _, (_, g) = Mooncake.value_and_gradient!!(cache, f, x0)
        @test g ≈ ForwardDiff.gradient(f, x0) rtol = 1e-10
    end
end
