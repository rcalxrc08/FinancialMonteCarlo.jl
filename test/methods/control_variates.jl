include("../testutils.jl")
using Random, Statistics

@testset "Control variates" begin
    model = BlackScholesProcess(sigma)
    Nstep_cv = 30
    variate_conf = MonteCarloConfiguration(100, Nstep_small)

    @testset "variate: $(nameof(typeof(variate)))" for variate in (Forward(T), AsianFloatingStrikeOption(T))
        mc_cv = MonteCarloConfiguration(Nsim_small, Nstep_cv; sampling = FinancialMonteCarlo.ControlVariates(variate, variate_conf))
        for payoff_ in (AsianFloatingStrikeOption(T), AsianFixedStrikeOption(T, K))
            @test isfinite(pricer(model, und, rfCurve, mc_cv, payoff_))
        end
    end

    @testset "reseed" begin
        eu = EuropeanOption(T, K)
        cv_method = FinancialMonteCarlo.ControlVariates(Forward(T), variate_conf)
        mc_reset = MonteCarloConfiguration(Nsim_small, Nstep_small; sampling = cv_method, seed = 1)
        price = pricer(model, und, rfCurve, mc_reset, eu)
        @test price == pricer(model, und, rfCurve, mc_reset, eu)
        mc_stream = MonteCarloConfiguration(Nsim_small, Nstep_small; sampling = cv_method, seed = 1, rng = MersenneTwister(1), reseed = false)
        @test pricer(model, und, rfCurve, mc_stream, eu) ≈ price
        @test pricer(model, und, rfCurve, mc_stream, eu) != price
    end

    @testset "the sampling of the main run keeps its checks" begin
        antithetic = FinancialMonteCarlo.ControlVariates(Forward(T), variate_conf, FinancialMonteCarlo.AntitheticMC())
        @test_throws ArgumentError MonteCarloConfiguration(1001, Nstep_small; sampling = antithetic)
        @test_throws ArgumentError MonteCarloConfiguration(1001, Nstep_small; sampling = FinancialMonteCarlo.AntitheticMC(), reseed = false)
        @test_throws ArgumentError MonteCarloConfiguration(1024, Nstep_small; sampling = FinancialMonteCarlo.ControlVariates(Forward(T), variate_conf, FinancialMonteCarlo.SobolMC()), mode = FinancialMonteCarlo.MultiThreading(2))
    end

    @testset "a Forward control variate reduces the variance" begin
        eu = EuropeanOption(T, K)
        seeds = 1:30
        plain = [pricer(model, und, rfCurve, MonteCarloConfiguration(2_000, 1; seed = s), eu) for s in seeds]
        cv = map(seeds) do s
            pilot = MonteCarloConfiguration(500, 1; seed = s)
            pricer(model, und, rfCurve, MonteCarloConfiguration(2_000, 1; sampling = FinancialMonteCarlo.ControlVariates(Forward(T), pilot), seed = s), eu)
        end
        @test std(cv) < std(plain) / 2
        @test abs(mean(cv) - BS_CALL_REF.eu) < 0.1
    end

    @testset "the exact price of the control (#75)" begin
        # an Asian option controlled by a European call: with the pilot mean as the price of the call the reduction mostly
        # cancels; with its closed form it is about halved
        asian, call = AsianFloatingStrikeOption(T), EuropeanOption(T, K)
        call_price = blsprice(S0, K, r, T, sigma, d)
        seeds = 1:40
        run(s; price = nothing) = pricer(model, und, rfCurve, MonteCarloConfiguration(2_000, 10; sampling = FinancialMonteCarlo.ControlVariates(call, MonteCarloConfiguration(2_000, 10; seed = 1000 + s); price = price), seed = s), asian)
        plain = [pricer(model, und, rfCurve, MonteCarloConfiguration(2_000, 10; seed = s), asian) for s in seeds]
        pilot, exact = run.(seeds), [run(s; price = call_price) for s in seeds]
        @test std(exact) < 0.75 * std(pilot)
        @test std(exact) < 0.75 * std(plain)
        @test abs(mean(exact) - mean(plain)) < 3 * std(plain) / sqrt(length(seeds))
        @test FinancialMonteCarlo.ControlVariates(call, MonteCarloConfiguration(10, 1); price = 1.0f0).price === 1.0f0
    end
end
