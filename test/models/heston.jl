include("../testutils.jl")

const HESTON_REF = (fwd = 98.8790717426047, eu = 7.9813557592924, barrier = 7.006564636309922, asian_float = 4.242256952013707)

@testset "Heston" begin
    sigma_zero, kappa, theta, lambda, rho_ = 0.2, 0.01, 0.03, 0.01, 0.0
    heston(σ, σ0) = HestonProcess(σ, σ0, lambda, kappa, rho_, theta, Underlying(S0, d))
    model = heston(sigma, sigma_zero)

    @testset "standard" begin
        check_prices(model, rfCurve, mc, standard_payoffs(), HESTON_REF; toll = 0.8)
    end

    @testset "antithetic" begin
        check_prices(model, rfCurve, mc_anti, standard_payoffs(), HESTON_REF; toll = 0.8)
    end

    @testset "antithetic agrees with standard for strong correlation" begin
        model_corr = HestonProcess(0.3, 0.2, 0.0, 2.0, -0.7, 0.04, Underlying(S0, d))
        mc_big = MonteCarloConfiguration(100_000, Nstep)
        mc_big_anti = MonteCarloConfiguration(100_000, Nstep, FinancialMonteCarlo.AntitheticMC())
        for payoff_ in (EuropeanOption(T, 80.0, false), EuropeanOption(T, 120.0))
            # about 7 standard errors of the standard estimator
            @test pricer(model_corr, rfCurve, mc_big_anti, payoff_) ≈ pricer(model_corr, rfCurve, mc_big, payoff_) atol = 0.1
        end
    end

    @testset "American with almost no path in the money at early dates" begin
        # fewer than 3 paths in the money make the regression singular: no early exercise there
        model_corr = HestonProcess(0.3, 0.2, 0.0, 2.0, -0.7, 0.04, Underlying(S0, d))
        @test isfinite(pricer(model_corr, rfCurve, MonteCarloConfiguration(5_000, 50), AmericanOption(T, 90.0, false)))
    end

    @testset "zero rate curve" begin
        rf_curve = ZeroRate([0.02; 0.02], T)
        @test abs(pricer(model, rf_curve, mc, Forward(T)) - 98.72567723404445) < 0.8
        @test abs(pricer(model, rf_curve, mc_anti, Forward(T)) - 98.72567723404445) < 0.8
    end

    @testset "parameter validation" begin
        @test_throws AssertionError heston(-sigma, sigma_zero)
        @test_throws AssertionError heston(sigma, -sigma_zero)
        @test_throws AssertionError HestonProcess(sigma, sigma_zero, lambda, kappa, -5.0, theta, Underlying(S0, d))
        @test_throws AssertionError HestonProcess(sigma, sigma_zero, 1e-16, 1e-16, rho_, theta, Underlying(S0, d))
    end
end
