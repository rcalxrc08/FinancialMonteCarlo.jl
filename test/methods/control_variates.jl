include("../testutils.jl")
using Random, Statistics

@testset "Control variates" begin
    model = BlackScholesProcess(sigma, Underlying(S0, d))
    Nstep_cv = 30
    variate_conf = MonteCarloConfiguration(100, Nstep_small)

    @testset "variate: $(nameof(typeof(variate)))" for variate in (Forward(T), AsianFloatingStrikeOption(T))
        mc_cv = MonteCarloConfiguration(Nsim_small, Nstep_cv; sampling = FinancialMonteCarlo.ControlVariates(variate, variate_conf))
        for payoff_ in (AsianFloatingStrikeOption(T), AsianFixedStrikeOption(T, K))
            @test isfinite(pricer(model, rfCurve, mc_cv, payoff_))
        end
    end

    @testset "reseed" begin
        eu = EuropeanOption(T, K)
        cv_method = FinancialMonteCarlo.ControlVariates(Forward(T), variate_conf)
        mc_reset = MonteCarloConfiguration(Nsim_small, Nstep_small; sampling = cv_method, seed = 1)
        price = pricer(model, rfCurve, mc_reset, eu)
        @test price == pricer(model, rfCurve, mc_reset, eu)
        mc_stream = MonteCarloConfiguration(Nsim_small, Nstep_small; sampling = cv_method, seed = 1, rng = MersenneTwister(1), reseed = false)
        @test pricer(model, rfCurve, mc_stream, eu) ≈ price
        @test pricer(model, rfCurve, mc_stream, eu) != price
    end

    @testset "the sampling of the main run keeps its checks" begin
        antithetic = FinancialMonteCarlo.ControlVariates(Forward(T), variate_conf, FinancialMonteCarlo.AntitheticMC())
        @test_throws AssertionError MonteCarloConfiguration(1001, Nstep_small; sampling = antithetic)
        @test_throws AssertionError MonteCarloConfiguration(1001, Nstep_small; sampling = FinancialMonteCarlo.AntitheticMC(), reseed = false)
        @test_throws ArgumentError MonteCarloConfiguration(1024, Nstep_small; sampling = FinancialMonteCarlo.ControlVariates(Forward(T), variate_conf, FinancialMonteCarlo.SobolMC()), mode = FinancialMonteCarlo.MultiThreading(2))
    end

    @testset "a Forward control variate reduces the variance" begin
        eu = EuropeanOption(T, K)
        seeds = 1:30
        plain = [pricer(model, rfCurve, MonteCarloConfiguration(2_000, 1; seed = s), eu) for s in seeds]
        cv = map(seeds) do s
            pilot = MonteCarloConfiguration(500, 1; seed = s)
            pricer(model, rfCurve, MonteCarloConfiguration(2_000, 1; sampling = FinancialMonteCarlo.ControlVariates(Forward(T), pilot), seed = s), eu)
        end
        @test std(cv) < std(plain) / 2
        @test abs(mean(cv) - BS_CALL_REF.eu) < 0.1
    end
end
