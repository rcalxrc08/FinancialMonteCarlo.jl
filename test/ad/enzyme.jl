include("../testutils.jl")
using Enzyme, ForwardDiff, Random

@testset "Enzyme" begin
    @testset "Kou European call: gradient matches ForwardDiff" begin
        # The default MersenneTwister: its scalar draws (jump times and sizes) are inactive in the extension.
        # x = [S0, σ, λ, p, λ₊, λ₋, r, d, T]
        x0 = [S0, sigma, 5.0, 0.3, 30.0, 20.0, r, d, T]
        mc_kou = MonteCarloConfiguration(100, Nstep_small)
        f(x) = pricer(KouProcess(x[2], x[3], x[4], x[5], x[6], Underlying(x[1], x[8])), ZeroRate(x[7]), mc_kou, EuropeanOption(x[9], K))
        g_ref = ForwardDiff.gradient(f, x0)
        g = Enzyme.gradient(Enzyme.Reverse, Enzyme.Const(f), x0)[1]
        @test g ≈ g_ref rtol = 1e-8
        @test collect(Enzyme.gradient(Enzyme.Forward, Enzyme.Const(f), x0)[1]) ≈ g_ref rtol = 1e-8
        @test 0 < g[1] < 1 # delta
        @test g[2] > 0 # vega
    end

    @testset "two-leg position built inside the function (#36, #40)" begin
        mc_xo = MonteCarloConfiguration(200, Nstep_small; rng = Xoshiro())
        f(x) = pricer(BlackScholesProcess(x[2], Underlying(x[1], d)), ZeroRate(r), mc_xo, 2 * EuropeanOption(T, x[3]) - EuropeanOption(T, x[3], false))
        x0 = [S0, sigma, K]
        @test Enzyme.gradient(Enzyme.Reverse, Enzyme.Const(f), x0)[1] ≈ ForwardDiff.gradient(f, x0) rtol = 1e-10
    end

    @testset "Black-Scholes, Xoshiro" begin
        x0 = [S0, sigma, r, d]
        mc_xo = MonteCarloConfiguration(2000, Nstep_small; rng = Xoshiro())
        f(x) = pricer(BlackScholesProcess(x[2], Underlying(x[1], x[4])), ZeroRate(x[3]), mc_xo, EuropeanOption(T, K))
        g = Enzyme.gradient(Enzyme.Reverse, Enzyme.Const(f), x0)[1]
        @test g ≈ ForwardDiff.gradient(f, x0) rtol = 1e-10
    end
end
