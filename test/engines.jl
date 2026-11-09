include("testutils.jl")
using Distributions, ForwardDiff, Statistics

# Engines are the building blocks of the models: they know nothing about rates or dividends.
@testset "Engines" begin
    mu_curve = ZeroRate([0.00, 0.02], T).r
    drift = 0.0

    @testset "BrownianMotion" begin
        @test_throws AssertionError BrownianMotion(-sigma, drift)
        @test_throws AssertionError BrownianMotion(-sigma, mu_curve)
        S = simulate(BrownianMotion(sigma, 0.01), mc_small_anti, T)
        @test size(S) == (Nsim_small, Nstep_small + 1)
        @test all(iszero, S[:, 1])
    end

    @testset "GeometricBrownianMotion" begin
        @test_throws AssertionError GeometricBrownianMotion(-sigma, drift, 1.0)
        @test_throws AssertionError GeometricBrownianMotion(-sigma, mu_curve, 1.0)
        # the paths are proportional to x0
        terminal(x0) = (FinancialMonteCarlo.set_seed!(mc_small); mean(simulate(GeometricBrownianMotion(sigma, 0.01, x0), mc_small, T)[:, end]))
        @test ForwardDiff.derivative(terminal, S0) ≈ terminal(S0) / S0
    end

    @testset "SubordinatedBrownianMotion" begin
        @test_throws AssertionError SubordinatedBrownianMotion(-sigma, drift, 0.01, InverseGaussian(1.0, 1.0))
        @test_throws AssertionError SubordinatedBrownianMotion(-sigma, mu_curve, 0.01, InverseGaussian(1.0, 1.0))
        S = simulate(SubordinatedBrownianMotion(sigma, 0.01, 0.1, Gamma(T / Nstep_small, 1.0)), mc_small, T)
        @test size(S) == (Nsim_small, Nstep_small + 1)
        @test all(iszero, S[:, 1])
    end
end

# The Longstaff-Schwartz regression solves its 3x3 normal equations on scalars, with the pivoting of the LAPACK LU.
@testset "Longstaff-Schwartz 3x3 solve" begin
    using Random
    solve(A, b) = FinancialMonteCarlo.solve_3x3(Tuple([A[1, :]; b[1]]), Tuple([A[2, :]; b[2]]), Tuple([A[3, :]; b[3]]))
    rng = MersenneTwister(1)
    for A in (rand(rng, 3, 3), [0.0 1 2; 1 0 3; 4 5 0], [1e-3 2 3; 4 5 6; 7 8 10], [1.0 2 3; 2 1 3; 3 3 1])
        b = rand(rng, 3)
        @test collect(solve(A, b)) ≈ A \ b rtol = 1e-12
        r(i) = (A[i, 1], A[i, 2], A[i, 3], b[i])
        @test @inferred(Union{Nothing, NTuple{3, Float64}}, FinancialMonteCarlo.solve_3x3(r(1), r(2), r(3))) isa NTuple{3, Float64}
    end
    @test isnothing(solve([1.0 2 3; 2 4 6; 1 1 1], ones(3))) # singular
    @test isnothing(solve(zeros(3, 3), ones(3)))
    # the normal equations of a regression on 1, S, S^2, with the types of a Dual price
    s = 90 .+ 20 .* rand(rng, 50)
    y = ForwardDiff.Dual.(rand(rng, 50), 1.0)
    A = [length(s) sum(s) sum(s .^ 2); sum(s) sum(s .^ 2) sum(s .^ 3); sum(s .^ 2) sum(s .^ 3) sum(s .^ 4)]
    x = FinancialMonteCarlo.solve_3x3((length(s) * one(s[1]), A[1, 2], A[1, 3], sum(y)), (A[2, 1], A[2, 2], A[2, 3], sum(y .* s)), (A[3, 1], A[3, 2], A[3, 3], sum(y .* s .^ 2)))
    @test ForwardDiff.value.(collect(x)) ≈ A \ ForwardDiff.value.([sum(y), sum(y .* s), sum(y .* s .^ 2)]) rtol = 1e-6
    @test ForwardDiff.partials.(collect(x), 1) ≈ A \ [length(s), sum(s), sum(s .^ 2)] rtol = 1e-6
    # through the pricer: with respect to the strike only the payoffs are dual numbers, not the paths
    model = BlackScholesProcess(sigma, Underlying(S0, d))
    dK = ForwardDiff.derivative(k -> pricer(model, ZeroRate(r), mc_small, AmericanOption(T, k, false)), K)
    @test 0 < dK < 1
end
