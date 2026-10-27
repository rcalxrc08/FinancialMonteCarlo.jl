include("../testutils.jl")
using DualNumbers

# Dual-valued parameters must flow through simulation and payoffs; the real part must match the plain prices.
@testset "DualNumbers" begin
    sigma_dual = dual(sigma, 1.0)

    @testset "Black-Scholes" begin
        model = BlackScholesProcess(sigma_dual, Underlying(S0, d))
        check_prices(model, rfCurve, mc, standard_payoffs(), BS_CALL_REF; toll = 0.8)
        check_prices(model, rfCurve, mc_anti, standard_payoffs(), BS_CALL_REF; toll = 0.6)
        @test pricer(model, rfCurve, mc, EuropeanOption(T, K)).epsilon > 0 # vega of a call
    end

    @testset "Heston" begin
        model = HestonProcess(sigma, dual(0.2, 1.0), 0.01, 0.01, 0.0, 0.03, Underlying(S0, d))
        heston_ref = (fwd = 98.8790717426047, eu = 7.9813557592924, barrier = 7.006564636309922, asian_float = 4.242256952013707)
        check_prices(model, rfCurve, mc, standard_payoffs(), heston_ref; toll = 0.8)
        check_prices(model, rfCurve, mc_anti, standard_payoffs(), heston_ref; toll = 0.8)
    end

    @testset "Kou" begin
        model = KouProcess(sigma_dual, 5.0, 0.3, 30.0, 20.0, Underlying(S0, d))
        kou_ref = (fwd = 99.41332633109904, eu = 10.347332240535199, barrier = 8.860123655599818, asian_float = 5.81798437145069)
        check_prices(model, rfCurve, mc, standard_payoffs(), kou_ref; toll = 0.8)
        check_prices(model, rfCurve, mc_anti, standard_payoffs(), kou_ref; toll = 0.8)
    end

    @testset "Variance Gamma" begin
        model = VarianceGammaProcess(sigma_dual, 0.01, 0.03, Underlying(S0, d))
        vg_ref = (fwd = 99.09861035102614, eu = 8.368903690692187, barrier = 7.469024475794258, asian_float = 4.779663836736272)
        check_prices(model, rfCurve, mc, standard_payoffs(), vg_ref; toll = 0.8)
        check_prices(model, rfCurve, mc_anti, standard_payoffs(), vg_ref; toll = 0.8)
        @test isfinite(pricer(model, rfCurve, MonteCarloConfiguration(Nsim, Nstep, FinancialMonteCarlo.SobolMode()), Forward(T)))
    end
end
