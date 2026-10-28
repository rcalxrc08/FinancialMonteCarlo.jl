include("../testutils.jl")
using Statistics

@testset "Control variates" begin
    model = BlackScholesProcess(sigma, Underlying(S0, d))
    Nstep_cv = 300
    variate_conf = MonteCarloConfiguration(1000, 100)

    @testset "variate: $(nameof(typeof(variate)))" for variate in (Forward(T), AsianFloatingStrikeOption(T))
        mc_cv = MonteCarloConfiguration(Nsim, Nstep_cv, FinancialMonteCarlo.ControlVariates(variate, variate_conf))
        for payoff_ in (AsianFloatingStrikeOption(T), AsianFixedStrikeOption(T, K))
            @test isfinite(pricer(model, rfCurve, mc_cv, payoff_))
        end
    end

    @testset "a Forward control variate reduces the variance" begin
        eu = EuropeanOption(T, K)
        seeds = 1:30
        plain = [pricer(model, rfCurve, MonteCarloConfiguration(2_000, 1, FinancialMonteCarlo.StandardMC(), FinancialMonteCarlo.SerialMode(s)), eu) for s in seeds]
        cv = map(seeds) do s
            pilot = MonteCarloConfiguration(500, 1, FinancialMonteCarlo.StandardMC(), FinancialMonteCarlo.SerialMode(s))
            pricer(model, rfCurve, MonteCarloConfiguration(2_000, 1, FinancialMonteCarlo.ControlVariates(Forward(T), pilot), FinancialMonteCarlo.SerialMode(s)), eu)
        end
        @test std(cv) < std(plain) / 2
        @test abs(mean(cv) - BS_CALL_REF.eu) < 0.1
    end
end
