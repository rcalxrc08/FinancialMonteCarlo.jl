include("../testutils.jl")

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
end
