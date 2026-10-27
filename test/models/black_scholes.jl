include("../testutils.jl")

@testset "Black-Scholes" begin
    model = BlackScholesProcess(sigma, Underlying(S0, d))

    @testset "parameters" begin
        m = BlackScholesProcess(0.1, Underlying(S0, d))
        @test FinancialMonteCarlo.get_parameters(m)[1] == 0.1
        m = FinancialMonteCarlo.set_parameters!(m, [sigma])
        @test m.σ == sigma
        @test FinancialMonteCarlo.get_parameters(Forward(T))[1] == T
        @test occursin("BlackScholesProcess", capture_stdout(() -> display(m)))
        @test_throws AssertionError BlackScholesProcess(-sigma, Underlying(S0, d))
    end

    @testset "calls" begin
        payoffs = merge(standard_payoffs(), (spot = Spot(), binary_eu = BinaryEuropeanOption(T, K)))
        check_prices(model, rfCurve, mc, payoffs, BS_CALL_REF; toll = 0.8)
    end

    @testset "calls (antithetic)" begin
        check_prices(model, rfCurve, mc_anti, standard_payoffs(), BS_CALL_REF; toll = 0.6)
    end

    @testset "puts" begin
        check_prices(model, rfCurve, mc, standard_put_payoffs(), BS_PUT_REF; toll = 0.6)
    end

    @testset "exotics" begin
        payoffs = (
            barrier_di = BarrierOptionDownIn(T, K, D),
            barrier_ui = BarrierOptionUpIn(T, K, D),
            barrier_uo = BarrierOptionUpOut(T, K, D),
            bermudan = BermudanOption(collect(0.2:0.1:T), K),
            binary_am = BinaryAmericanOption(T, K),
            binary_eu = BinaryEuropeanOption(T, K),
            double_barrier = DoubleBarrierOption(T, K, K / 10.0, 1.2 * K),
        )
        check_prices(model, rfCurve, mc, payoffs)
    end

    @testset "smiles" begin
        K_vec = [80.0, 90.0, 100.0, 110.0, 120.0]
        isCall_vec = [true, true, true, false, true]
        for smile in (FinancialMonteCarlo.EuropeanOptionSmile(T, K_vec, isCall_vec), FinancialMonteCarlo.BinaryEuropeanOptionSmile(T, K_vec, isCall_vec))
            prices = pricer(model, rfCurve, mc, smile)
            @test length(prices) == length(K_vec)
            @test all(isfinite, prices)
        end
    end

    @testset "time dependent rates and dividends" begin
        mc_large = MonteCarloConfiguration(100_000, Nstep)
        model_curve = BlackScholesProcess(sigma, Underlying(S0, FinancialMonteCarlo.Curve([0.019, 0.019, 0.02], T)))
        rf_curve = ZeroRate([0.00, 0.02], T)
        @test isfinite(pricer(model_curve, rf_curve, MonteCarloConfiguration(100_000, Nstep, FinancialMonteCarlo.AntitheticMC()), Forward(T)))
        payoffs = merge(standard_payoffs(), (binary_eu = BinaryEuropeanOption(T, K),))
        check_prices(model_curve, rf_curve, mc_large, payoffs, BS_CALL_REF; toll = 1.8)
    end
end
