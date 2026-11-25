include("../testutils.jl")
using DualNumbers

# Finite-difference greeks are checked against dual-number derivatives on the same seed.
@testset "Greeks" begin
    @testset "delta" begin
        name_1, name_2 = "ENJ", "ABPL"
        basket_name = (name_1, name_2, "TESL")
        rho_1 = [1.0 0.0 0.1; 0.0 1.0 0.0; 0.1 0.0 1.0]
        model_enj = BlackScholesProcess(sigma)
        und_dual = Underlying(dual(S0, 1.0), d)
        basket = GaussianCopulaNVariateProcess(rho_1, model_enj, model_enj, model_enj)
        am = AmericanOption(T, K)
        portfolio_basket = basket_name → EuropeanOptionND(T, K)
        portfolio_single = name_1 → am
        # the spot of name_1 is a dual number: the derivative of the price in it is the delta
        basket_market(u1) = basket_name → (basket, (u1, und, und))

        @test occursin("GaussianCopulaNVariateProcess", capture_stdout(() -> display(basket)))

        delta_basket = FinancialMonteCarlo.delta(basket_market(und), rfCurve, mc, portfolio_basket, name_1, 1e-7)
        delta_single = FinancialMonteCarlo.delta(name_1 → (model_enj, und), rfCurve, mc, portfolio_single, name_1, 1e-7)
        dual_basket = pricer(basket_market(und_dual), rfCurve, mc, portfolio_basket).epsilon
        dual_single = pricer(name_1 → (model_enj, und_dual), rfCurve, mc, portfolio_single).epsilon

        @test abs(dual_basket - delta_basket) < 1e-4
        @test abs(dual_single - delta_single) < 1e-4
        @test abs(dual_single - FinancialMonteCarlo.delta(model_enj, und, rfCurve, mc, am, 1e-7)) < 1e-4
        @test abs(dual_single - FinancialMonteCarlo.delta(model_enj, und, rfCurve, mc, am * 1.0, 1e-7)) < 1e-4
        @test abs(dual_single - FinancialMonteCarlo.delta(model_enj, und, rfCurve, mc, [am], 1e-7)[1]) < 1e-4

        @test_throws ArgumentError FinancialMonteCarlo.delta(name_1 → (model_enj, und), rfCurve, mc, portfolio_single, "ENJ_ABPL", 1e-7)
        @test FinancialMonteCarlo.delta(name_1 → (model_enj, und), rfCurve, mc, portfolio_single, "MISSING", 1e-7) == 0.0
    end

    @testset "rho" begin
        dr = 1e-7
        der_toll = 1e-2
        model = BlackScholesProcess(sigma)
        rf_dual = ZeroRate(dual(r, 1.0))

        @testset "$(name)" for (name, payoff_, cfg, toll) in [("european", EuropeanOption(T, K), mc, 0.8), ("american", AmericanOption(T, K), mc, 0.8), ("european antithetic", EuropeanOption(T, K), mc_anti, 0.6)]
            price = pricer(model, und, rf_dual, cfg, payoff_)
            @test abs(price - (payoff_ isa AmericanOption ? BS_CALL_REF.am : BS_CALL_REF.eu)) < toll
            @test abs(price.epsilon - rho(model, und, rfCurve, cfg, payoff_, dr)) < der_toll
            @test abs(price.epsilon - rho(model, und, rfCurve, cfg, [payoff_], dr)[1]) < der_toll
        end

        @testset "parallel shift of a curve" begin
            flat_curve = ZeroRate(FinancialMonteCarlo.Curve([r, r, r], T))
            @test rho(model, und, flat_curve, mc, EuropeanOption(T, K), 1e-5) ≈ rho(model, und, rfCurve, mc, EuropeanOption(T, K), 1e-5) rtol = 1e-3
        end
    end
end
