include("testutils.jl")
using Distributions

# Engines are the building blocks of the models: they know nothing about rates or dividends.
@testset "Engines" begin
    mu_curve = ZeroRate([0.00, 0.02], T).r
    drift = 0.0

    @testset "BrownianMotion" begin
        @test_throws AssertionError BrownianMotion(-sigma, drift)
        @test_throws AssertionError BrownianMotion(-sigma, mu_curve)
        S = simulate(BrownianMotion(sigma, 0.01), mc_anti, T)
        @test size(S) == (Nsim, Nstep + 1)
        @test all(iszero, S[:, 1])
    end

    @testset "GeometricBrownianMotion" begin
        @test_throws AssertionError GeometricBrownianMotion(-sigma, drift, 1.0)
        @test_throws AssertionError GeometricBrownianMotion(-sigma, mu_curve, 1.0)
    end

    @testset "SubordinatedBrownianMotion" begin
        @test_throws AssertionError SubordinatedBrownianMotion(-sigma, drift, 0.01, InverseGaussian(1.0, 1.0))
        @test_throws AssertionError SubordinatedBrownianMotion(-sigma, mu_curve, 0.01, InverseGaussian(1.0, 1.0))
    end
end
