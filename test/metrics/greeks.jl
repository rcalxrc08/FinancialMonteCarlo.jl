include("../testutils.jl")
using DualNumbers

# Finite-difference greeks are checked against dual-number derivatives on the same seed.
@testset "Greeks" begin
    @testset "delta" begin
        name_1, name_2 = "ENJ", "ABPL"
        basket_name = name_1 * "_" * name_2 * "_TESL"
        rho_1 = [1.0 0.0 0.1; 0.0 1.0 0.0; 0.1 0.0 1.0]
        bs(S) = BlackScholesProcess(sigma, Underlying(S, d))
        model_enj, model_abpl, model_tesl = bs(S0), bs(S0), bs(S0)
        model_enj_dual = bs(dual(S0, 1.0))
        basket = GaussianCopulaNVariateProcess(rho_1, model_enj, model_abpl, model_tesl)
        basket_dual = GaussianCopulaNVariateProcess(rho_1, model_enj_dual, model_abpl, model_tesl)
        am = AmericanOption(T, K)
        portfolio_basket = basket_name → EuropeanOptionND(T, K)
        portfolio_single = name_1 → am

        @test occursin("GaussianCopulaNVariateProcess", capture_stdout(() -> display(basket)))

        delta_basket = FinancialMonteCarlo.delta(basket_name → basket, rfCurve, mc, portfolio_basket, name_1, 1e-7)
        delta_single = FinancialMonteCarlo.delta(name_1 → model_enj, rfCurve, mc, portfolio_single, name_1, 1e-7)
        dual_basket = pricer(basket_name → basket_dual, rfCurve, mc, portfolio_basket).epsilon
        dual_single = pricer(name_1 → model_enj_dual, rfCurve, mc, portfolio_single).epsilon

        @test abs(dual_basket - delta_basket) < 1e-4
        @test abs(dual_single - delta_single) < 1e-4
        @test abs(dual_single - FinancialMonteCarlo.delta(model_enj, rfCurve, mc, am, 1e-7)) < 1e-4
        @test abs(dual_single - FinancialMonteCarlo.delta(model_enj, rfCurve, mc, am * 1.0, 1e-7)) < 1e-4
        @test abs(dual_single - FinancialMonteCarlo.delta(model_enj, rfCurve, mc, [am], 1e-7)[1]) < 1e-4

        @test_throws AssertionError FinancialMonteCarlo.delta(name_1 → model_enj, rfCurve, mc, portfolio_single, basket_name, 1e-7)
        @test FinancialMonteCarlo.delta(name_1 → model_enj, rfCurve, mc, portfolio_single, "MISSING", 1e-7) == 0.0
    end

    @testset "rho" begin
        dr = 1e-7
        der_toll = 1e-2
        model = BlackScholesProcess(sigma, Underlying(S0, d))
        rf_dual = ZeroRate(dual(r, 1.0))

        @testset "$(name)" for (name, payoff_, cfg, toll) in [("european", EuropeanOption(T, K), mc, 0.8), ("american", AmericanOption(T, K), mc, 0.8), ("european antithetic", EuropeanOption(T, K), mc_anti, 0.6)]
            price = pricer(model, rf_dual, cfg, payoff_)
            @test abs(price - (payoff_ isa AmericanOption ? BS_CALL_REF.am : BS_CALL_REF.eu)) < toll
            @test abs(price.epsilon - rho(model, rfCurve, cfg, payoff_, dr)) < der_toll
            @test abs(price.epsilon - rho(model, rfCurve, cfg, [payoff_], dr)[1]) < der_toll
        end
    end
end
