include("../testutils.jl")
using FinancialToolbox: blsprice

@testset "ShiftedProcess" begin
    und = Underlying(S0, d)
    α = 0.3
    model = ShiftedProcess(BlackScholesProcess(sigma, und), α)

    @testset "construction and parameters" begin
        @test model.underlying === und
        @test FinancialMonteCarlo.get_parameters(model) == [sigma, α]
        new_model = FinancialMonteCarlo.set_parameters!(model, [0.3, 0.1])
        @test (new_model.model.σ, new_model.α) == (0.3, 0.1)
        @test_throws DimensionMismatch FinancialMonteCarlo.set_parameters!(model, [0.3])
        @test_throws DomainError ShiftedProcess(BlackScholesProcess(sigma, und), 1.0)
        @test occursin("ShiftedProcess", capture_stdout(() -> display(model)))
    end

    @testset "displaced Black-Scholes against the formula" begin
        # S_T - α F is a GBM started from (1 - α) S0: a call on S_T is a call on it with strike K - α F
        F = S0 * exp((r - d) * T)
        mc_sobol = MonteCarloConfiguration(2^17, 1; sampling = FinancialMonteCarlo.SobolMC())
        for strike in (80.0, K, 120.0), isCall in (true, false)
            @test abs(pricer(model, rfCurve, mc_sobol, EuropeanOption(T, strike, isCall)) - blsprice((1 - α) * S0, strike - α * F, r, T, sigma, d, isCall)) < 1e-2
        end
        refs = (fwd = S0 * exp(-d * T), eu = blsprice((1 - α) * S0, K - α * F, r, T, sigma, d))
        check_prices(model, rfCurve, mc, standard_payoffs(), refs; toll = 0.6)
        check_prices(model, rfCurve, mc_anti, standard_payoffs(), refs; toll = 0.6)
    end

    @testset "no shift is the model itself" begin
        @test pricer(ShiftedProcess(BlackScholesProcess(sigma, und), 0.0), rfCurve, mc_small, EuropeanOption(T, K)) ≈ pricer(BlackScholesProcess(sigma, und), rfCurve, mc_small, EuropeanOption(T, K))
    end

    @testset "a shifted Bachelier is the Bachelier" begin
        # (1 - α) S0 + σ W + α S0 = S0 + σ W, path by path
        bachelier = BachelierProcess(20.0, und)
        @test pricer(ShiftedProcess(bachelier, -0.2), rfCurve, mc_small, AsianFixedStrikeOption(T, K)) ≈ pricer(bachelier, rfCurve, mc_small, AsianFixedStrikeOption(T, K))
    end

    @testset "$(nameof(typeof(inner)))" for inner in (HestonProcess(0.5, 0.2, 0.0, 2.0, -0.7, 0.04, und), KouProcess(0.2, 3.0, 0.3, 9.0, 5.0, und), MertonProcess(0.2, 3.0, -0.05, 0.1, und), VarianceGammaProcess(0.2, -0.1, 0.2, und), NormalInverseGaussianProcess(0.2, -0.1, 0.2, und))
        shifted = ShiftedProcess(inner, 0.2)
        # the forward is unchanged, and on the same seed the shifted call is cheaper: the payoff is convex
        @test abs(pricer(shifted, rfCurve, mc_anti, Forward(T)) - S0 * exp(-d * T)) < 0.5
        @test pricer(shifted, rfCurve, mc, EuropeanOption(T, K)) < pricer(inner, rfCurve, mc, EuropeanOption(T, K))
        for config in (mc_small, mc_small_anti, MonteCarloConfiguration(Nsim_small, Nstep_small; sampling = FinancialMonteCarlo.SobolMC()))
            check_prices(shifted, rfCurve, config, standard_payoffs())
        end
        check_prices(shifted, ZeroRate([0.01, 0.02], T), mc_small, standard_payoffs())
    end

    @testset "Float32" begin
        model32 = ShiftedProcess(BlackScholesProcess(0.2f0, Underlying(100.0f0, 0.01f0)), 0.3f0)
        @test pricer(model32, ZeroRate(0.02f0), mc_small, EuropeanOption(1.0f0, 100.0f0)) isa Float32
    end
end
