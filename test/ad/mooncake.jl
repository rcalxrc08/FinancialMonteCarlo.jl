include("../testutils.jl")
using Mooncake, ForwardDiff, Random

@testset "Mooncake" begin
    mc_xo = MonteCarloConfiguration(100, Nstep_small; rng = Xoshiro())

    @testset "Kou European call: gradient matches ForwardDiff" begin
        # x = [S0, σ, λ, p, λ₊, λ₋, r, d, T]
        x0 = [S0, sigma, 5.0, 0.3, 30.0, 20.0, r, d, T]
        f(x) = pricer(KouProcess(x[2], x[3], x[4], x[5], x[6]), Underlying(x[1], x[8]), ZeroRate(x[7]), mc_xo, EuropeanOption(x[9], K))
        cache = Mooncake.prepare_gradient_cache(f, x0)
        price, (_, g) = Mooncake.value_and_gradient!!(cache, f, x0)
        @test price ≈ f(x0)
        @test g ≈ ForwardDiff.gradient(f, x0) rtol = 1e-8
        @test 0 < g[1] < 1 # delta
        @test g[2] > 0 # vega
    end

    @testset "two-leg position built inside the function (#36)" begin
        f(x) = pricer(BlackScholesProcess(x[2]), Underlying(x[1], d), ZeroRate(r), mc_xo, 2 * EuropeanOption(T, x[3]) - EuropeanOption(T, x[3], false))
        x0 = [S0, sigma, K]
        cache = Mooncake.prepare_gradient_cache(f, x0)
        _, (_, g) = Mooncake.value_and_gradient!!(cache, f, x0)
        @test g ≈ ForwardDiff.gradient(f, x0) rtol = 1e-8
    end

    @testset "the rules of the kernels: reverse and forward mode match ForwardDiff" begin
        # every kernel with a rule in the extension (bm_step!, exp_paths!, affine_carry!, add_carry!, heston_step!,
        # lnm_step!, add_jump_sizes!, path_range), antithetic, with a rate curve and three payoffs on the same paths
        x0 = [S0, sigma, r, d]
        mc_anti_xo = MonteCarloConfiguration(200, Nstep_small; sampling = FinancialMonteCarlo.AntitheticMC(), rng = Xoshiro())
        legs = EuropeanOption(T, K) + AsianFixedStrikeOption(T, K) + BarrierOptionDownOut(T, K, 85.0)
        models = (x -> BlackScholesProcess(x[2]), x -> BachelierProcess(100 * x[2]), x -> HestonProcess(0.5, x[2], 0.0, 2.0, -0.7, 0.04), x -> KouProcess(x[2], 3.0, 0.3, 30.0, 20.0), x -> MertonProcess(x[2], 3.0, -0.05, 0.1), x -> LogNormalMixture((x[2], 0.3), (0.5,)), x -> ShiftedProcess(BlackScholesProcess(x[2]), 0.2))
        @testset "model $k" for (k, build) in enumerate(models)
            f(x) = pricer(build(x), Underlying(x[1], x[4]), ZeroRate(FinancialMonteCarlo.Curve([x[3], 2 * x[3]], T)), mc_anti_xo, legs)
            g_ref = ForwardDiff.gradient(f, x0)
            cache = Mooncake.prepare_gradient_cache(f, x0)
            _, (_, g) = Mooncake.value_and_gradient!!(cache, f, x0)
            @test g ≈ g_ref rtol = 1e-10
            forward = Mooncake.prepare_derivative_cache(f, x0)
            @test [Mooncake.value_and_derivative!!(forward, (f, Mooncake.zero_tangent(f)), (x0, Float64.(1:4 .== i)))[2] for i = 1:4] ≈ g_ref rtol = 1e-10
        end
    end

    @testset "default rng: seeding has a rule, refills go through" begin
        # 2000 paths draw past the 1002 numbers cached by the MersenneTwister.
        x0 = [S0, sigma, r, d]
        mc_default = MonteCarloConfiguration(2000, Nstep_small)
        f(x) = pricer(BlackScholesProcess(x[2]), Underlying(x[1], x[4]), ZeroRate(x[3]), mc_default, EuropeanOption(T, K))
        cache = Mooncake.prepare_gradient_cache(f, x0)
        _, (_, g) = Mooncake.value_and_gradient!!(cache, f, x0)
        @test g ≈ ForwardDiff.gradient(f, x0) rtol = 1e-10
    end
end
