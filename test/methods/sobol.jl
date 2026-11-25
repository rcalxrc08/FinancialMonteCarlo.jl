include("../testutils.jl")

@testset "Sobol quasi Monte Carlo" begin
    mc_sobol = MonteCarloConfiguration(Nsim, Nstep; sampling = FinancialMonteCarlo.SobolMC())

    # sobol_skip! sets the state in closed form: the same state, and points, as Sobol.skip drawing every skipped point
    @testset "skip in closed form, dimension $dim" for dim in (1, 2, 7, 100)
        Sobol = FinancialMonteCarlo.Sobol
        for n in (0, 1, 2, 3, 8, 1023, 1024, 409_600)
            reference, fast = Sobol.SobolSeq(dim), Sobol.SobolSeq(dim)
            skip(reference, n)
            FinancialMonteCarlo.sobol_skip!(fast, n)
            @test (fast.x, fast.b, fast.n) == (reference.x, reference.b, reference.n)
            @test all(Sobol.next!(fast) == Sobol.next!(reference) for _ = 1:50)
        end
    end

    @testset "Black-Scholes" begin
        model = BlackScholesProcess(sigma)
        payoffs = merge(standard_payoffs(), (binary_eu = BinaryEuropeanOption(T, K),))
        check_prices(model, und, rfCurve, mc_sobol, payoffs, BS_CALL_REF; toll = 0.8)
        model_curve = BlackScholesProcess(sigma)
        @test isfinite(pricer(model_curve, Underlying(S0, FinancialMonteCarlo.Curve([0.009999, 0.01], T)), rfCurve, mc_sobol, Forward(T)))
    end

    @testset "Heston" begin
        model = HestonProcess(sigma, 0.2, 0.01, 0.01, 0.0, 0.03)
        heston_ref = (fwd = 98.8790717426047, eu = 7.9813557592924, barrier = 7.006564636309922, asian_float = 4.242256952013707)
        check_prices(model, und, rfCurve, mc_sobol, standard_payoffs(), heston_ref; toll = 0.8)
        @test isfinite(pricer(model, und, ZeroRate(FinancialMonteCarlo.ImpliedCurve([0.01], [0.5])), mc_sobol, Forward(T)))
        # a flat rate curve gives the same forward as the scalar rate
        @test abs(pricer(model, und, ZeroRate(FinancialMonteCarlo.Curve([r, r, r], T)), mc_sobol, Forward(T)) - S0 * exp(-d * T)) < 0.2
    end

    @testset "Variance Gamma" begin
        model = VarianceGammaProcess(sigma, 0.01, 0.03)
        vg_ref = (fwd = 99.09861035102614, eu = 8.368903690692187, barrier = 7.469024475794258, asian_float = 4.779663836736272)
        check_prices(model, und, rfCurve, mc_sobol, standard_payoffs(), vg_ref; toll = 0.8)
        model_curve = VarianceGammaProcess(sigma, 0.01, 0.03)
        @test abs(pricer(model_curve, Underlying(S0, FinancialMonteCarlo.Curve([0.009999, 0.01], T)), rfCurve, mc_sobol, Forward(T)) - vg_ref.fwd) < 0.8
    end
end
