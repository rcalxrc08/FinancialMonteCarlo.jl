include("../testutils.jl")
using Distributions, Statistics

const MERTON_REF = (fwd = 99.1188767166039, eu = 9.084327245917533, barrier = 7.880881290426765, asian_float = 5.129020349580892)

@testset "Merton" begin
    lam, mu1, sigma1 = 5.0, 0.03, 0.02
    model = MertonProcess(sigma, lam, mu1, sigma1, Underlying(S0, d))

    @testset "standard" begin
        check_prices(model, rfCurve, mc, standard_payoffs(), MERTON_REF; toll = 0.8)
    end

    @testset "antithetic" begin
        check_prices(model, rfCurve, mc_anti, standard_payoffs(), MERTON_REF; toll = 0.6)
    end

    @testset "sobol" begin
        mc_sobol = MonteCarloConfiguration(Nsim, Nstep; sampling = FinancialMonteCarlo.SobolMC())
        check_prices(model, rfCurve, mc_sobol, standard_payoffs(), MERTON_REF; toll = 0.6)
        # the jumps are quasi random too: the forward is close to its exact value
        @test abs(pricer(model, rfCurve, mc_sobol, Forward(T)) - S0 * exp(-d * T)) < 0.05
        # a path with more jumps than the Sobol point has room for draws them from the rng
        big = MertonProcess(sigma, 2000.0, 0.0, 0.001, Underlying(S0, d))
        @test isfinite(pricer(big, rfCurve, MonteCarloConfiguration(64, 2; sampling = FinancialMonteCarlo.SobolMC()), Forward(T)))
    end

    @testset "antithetic jumps" begin
        # the jumps of the path i + Nsim / 2 mirror those of the path i: with no diffusion to speak of, the log
        # returns of the two halves are the jump sums from u and 1 - u
        model_jumps = MertonProcess(1e-8, lam, 0.0, sigma1, Underlying(S0, 0.0))
        S = simulate(model_jumps, ZeroRate(0.0), mc_small_anti, T)
        half = div(Nsim_small, 2)
        jumps = log.(S[:, end] ./ S0)
        # symmetric jump law and counts N(u), N(1 - u): the pairs are negatively correlated
        @test cor(jumps[1:half], jumps[(half+1):end]) < 0
    end

    @testset "jump counts by inversion" begin
        for λT in (1e-3, 0.5, 5.0, 50.0), u in (0.0, 1e-6, 0.3, 0.5, 0.9, 1 - 1e-6, 1 - 1e-12)
            J = FinancialMonteCarlo.JumpUniforms(λT, Float64)
            @test FinancialMonteCarlo.jump_count(J, u) == quantile(Poisson(λT), u)
        end
        J = FinancialMonteCarlo.JumpUniforms(5.0, Float64)
        @test J.M == quantile(Poisson(5.0), 1 - sqrt(eps()))
        @test FinancialMonteCarlo.JumpUniforms(5.0f0, Float32).M < J.M
        @test FinancialMonteCarlo.JumpUniforms(2000.0, Float64).M > 2000
        @test FinancialMonteCarlo.jump_quantile(model, 0.3) ≈ mu1 + sigma1 * quantile(Normal(), 0.3)
    end

    @testset "parameter validation" begin
        @test_throws AssertionError MertonProcess(-sigma, lam, mu1, sigma1, Underlying(S0, d))
        @test_throws AssertionError MertonProcess(sigma, lam, mu1, -sigma1, Underlying(S0, d))
        @test_throws AssertionError MertonProcess(sigma, -lam, mu1, sigma1, Underlying(S0, d))
    end
end
