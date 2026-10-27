include("../testutils.jl")

@testset "Sobol quasi Monte Carlo" begin
    mc_sobol = MonteCarloConfiguration(Nsim, Nstep, FinancialMonteCarlo.SobolMode())

    @testset "Black-Scholes" begin
        model = BlackScholesProcess(sigma, Underlying(S0, d))
        payoffs = merge(standard_payoffs(), (binary_eu = BinaryEuropeanOption(T, K),))
        check_prices(model, rfCurve, mc_sobol, payoffs, BS_CALL_REF; toll = 0.8)
        model_curve = BlackScholesProcess(sigma, Underlying(S0, FinancialMonteCarlo.Curve([0.009999, 0.01], T)))
        @test isfinite(pricer(model_curve, rfCurve, mc_sobol, Forward(T)))
    end

    @testset "Heston" begin
        model = HestonProcess(sigma, 0.2, 0.01, 0.01, 0.0, 0.03, Underlying(S0, d))
        heston_ref = (fwd = 98.8790717426047, eu = 7.9813557592924, barrier = 7.006564636309922, asian_float = 4.242256952013707)
        check_prices(model, rfCurve, mc_sobol, standard_payoffs(), heston_ref; toll = 0.8)
        @test isfinite(pricer(model, FinancialMonteCarlo.ImpliedZeroRate([0.01], [0.5]), mc_sobol, Forward(T)))
    end

    @testset "Variance Gamma" begin
        model = VarianceGammaProcess(sigma, 0.01, 0.03, Underlying(S0, d))
        vg_ref = (fwd = 99.09861035102614, eu = 8.368903690692187, barrier = 7.469024475794258, asian_float = 4.779663836736272)
        check_prices(model, rfCurve, mc_sobol, standard_payoffs(), vg_ref; toll = 0.8)
        model_curve = VarianceGammaProcess(sigma, 0.01, 0.03, Underlying(S0, FinancialMonteCarlo.Curve([0.009999, 0.01], T)))
        @test abs(pricer(model_curve, rfCurve, mc_sobol, Forward(T)) - vg_ref.fwd) < 0.8
    end
end
