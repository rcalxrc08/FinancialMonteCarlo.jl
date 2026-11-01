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
        check_prices(model, rfCurve, mc_small, payoffs)
    end

    @testset "Bermudan and American engines" begin
        put_eu, put_am = EuropeanOption(T, K, false), AmericanOption(T, K, false)
        @test pricer(model, rfCurve, mc_small, BermudanOption([T], K, false)) ≈ pricer(model, rfCurve, mc_small, put_eu)
        # exercisable at every step of the grid, a Bermudan is the American on the same paths
        every_step = collect(T / Nstep_small:T / Nstep_small:T)
        @test pricer(model, rfCurve, mc_small, BermudanOption(every_step, K, false)) ≈ pricer(model, rfCurve, mc_small, put_am)
        @test abs(pricer(model, rfCurve, mc, put_am) - BS_PUT_REF.am) < 0.3
    end

    @testset "smiles" begin
        K_vec = [80.0, 90.0, 100.0, 110.0, 120.0]
        isCall_vec = [true, true, true, false, true]
        for smile in (FinancialMonteCarlo.EuropeanOptionSmile(T, K_vec, isCall_vec), FinancialMonteCarlo.BinaryEuropeanOptionSmile(T, K_vec, isCall_vec))
            options = FinancialMonteCarlo.smile_to_options(smile)
            @test pricer(model, rfCurve, mc_small, smile) ≈ pricer(model, rfCurve, mc_small, options)
            @test variance(model, rfCurve, mc_small, smile) ≈ variance(model, rfCurve, mc_small, options)
            @test delta(model, rfCurve, mc_small, smile) ≈ delta(model, rfCurve, mc_small, options)
            @test rho(model, rfCurve, mc_small, smile) ≈ rho(model, rfCurve, mc_small, options)
            @test all(collect.(confinter(model, rfCurve, mc_small, smile)) .≈ collect.(confinter(model, rfCurve, mc_small, options)))

            @test FinancialMonteCarlo.payout(95.0, smile) == [FinancialMonteCarlo.payout(95.0, option) for option in options]

            S = simulate(model, rfCurve, mc_small, 2 * T)
            payoffs = payoff(S, smile, rfCurve, mc_small, 2 * T)
            @test size(payoffs) == (Nsim_small, length(K_vec))
            @test all(payoffs[:, i] ≈ payoff(S, options[i], rfCurve, mc_small, 2 * T) for i in eachindex(options))
        end

        batches = MonteCarloConfiguration(Nsim_small, Nstep_small, FinancialMonteCarlo.StandardMC(), FinancialMonteCarlo.MultiThreading(FinancialMonteCarlo.SerialMode(), 1, [1]))
        prices = pricer(model, rfCurve, batches, FinancialMonteCarlo.EuropeanOptionSmile(T, K_vec, isCall_vec))
        @test prices isa Vector{Float64}
        @test length(prices) == length(K_vec)
    end

    @testset "time dependent rates and dividends" begin
        model_curve = BlackScholesProcess(sigma, Underlying(S0, FinancialMonteCarlo.Curve([0.019, 0.019, 0.02], T)))
        rf_curve = ZeroRate([0.00, 0.02], T)
        @test isfinite(pricer(model_curve, rf_curve, mc_small_anti, Forward(T)))
        payoffs = merge(standard_payoffs(), (binary_eu = BinaryEuropeanOption(T, K),))
        check_prices(model_curve, rf_curve, mc, payoffs, BS_CALL_REF; toll = 1.8)
    end
end
