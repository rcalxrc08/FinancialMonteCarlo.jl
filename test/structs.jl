include("testutils.jl")

# A payoff with no `payoff` method, to check the generic fallback errors out.
struct UnpricedPayoff <: FinancialMonteCarlo.AbstractPayoff{Float64}
    T::Float64
end

@testset "Structs" begin
    Kneg = -K
    Tneg = -T

    @testset "MonteCarloConfiguration" begin
        @test MonteCarloConfiguration(Nsim, Nstep, 1).parallelMode.seed == 1
        @test MonteCarloConfiguration(Nsim, Nstep, FinancialMonteCarlo.SerialMode(1)).parallelMode.seed == 1
        @test_throws AssertionError MonteCarloConfiguration(-Nsim, Nstep)
        @test_throws AssertionError MonteCarloConfiguration(Nsim, -Nstep)
        @test_throws AssertionError MonteCarloConfiguration(Nsim + 1, Nstep, FinancialMonteCarlo.AntitheticMC())
    end

    @testset "Underlying and curves" begin
        rf_curve = ZeroRate([0.00, 0.02], T)
        @test FinancialMonteCarlo.Curve([0.0, 0.1], [0.0, 0.1]) isa FinancialMonteCarlo.CurveType
        @test_throws AssertionError Underlying(-S0)
        @test_throws AssertionError Underlying(-S0, rf_curve.r)

        # trapezoidal inside the nodes, flat outside them
        c = ZeroRate([0.01, 0.02, 0.05], 1.0).r # nodes at 0, 0.5, 1
        @test FinancialMonteCarlo.integral(c, 1.0) ≈ 0.025
        @test FinancialMonteCarlo.integral(c, 2.0) ≈ 0.025 + 0.05
        c_late = FinancialMonteCarlo.Curve([0.03, 0.05], [0.5, 1.0])
        @test FinancialMonteCarlo.integral(c_late, 0.25) ≈ 0.03 * 0.25
        @test FinancialMonteCarlo.integral(c_late, 1.0) ≈ 0.03 * 0.5 + 0.04 * 0.5
        @test FinancialMonteCarlo.integral(FinancialMonteCarlo.Curve([0.03], [1.0]), 2.0) ≈ 0.06
    end

    @testset "market data and portfolios" begin
        model = GaussianCopulaNVariateProcess(BlackScholesProcess(sigma, Underlying(S0, d)), BlackScholesProcess(sigma, Underlying(S0, d)), 0.0)
        mkt = "eni_aap" → model
        @test_throws AssertionError mkt + mkt
        @test_throws AssertionError "ENJ" → model # multivariate key must be INDEX1_INDEX2
        @test_throws AssertionError "c_i" → Forward(1.0) # single name key cannot contain '_'
    end

    @testset "payoffs: negative maturity" begin
        @test_throws AssertionError Forward(Tneg)
        @test_throws AssertionError EuropeanOption(Tneg, K)
        @test_throws AssertionError EuropeanOptionND(Tneg, K)
        @test_throws AssertionError AmericanOption(Tneg, K)
        @test_throws AssertionError BermudanOption([Tneg], K)
        @test_throws AssertionError BinaryAmericanOption(Tneg, K)
        @test_throws AssertionError BinaryEuropeanOption(Tneg, K)
        @test_throws AssertionError BarrierOptionDownOut(Tneg, K, D)
        @test_throws AssertionError BarrierOptionUpOut(Tneg, K, D)
        @test_throws AssertionError BarrierOptionUpIn(Tneg, K, D)
        @test_throws AssertionError BarrierOptionDownIn(Tneg, K, D)
        @test_throws AssertionError AsianFloatingStrikeOption(Tneg)
        @test_throws AssertionError AsianFixedStrikeOption(Tneg, K)
        @test_throws AssertionError DoubleBarrierOption(Tneg, K, K * 10, D)
    end

    @testset "payoffs: negative strike or barrier" begin
        @test_throws AssertionError EuropeanOption(T, Kneg)
        @test_throws AssertionError EuropeanOptionND(T, Kneg)
        @test_throws AssertionError AmericanOption(T, Kneg)
        @test_throws AssertionError BermudanOption([T], Kneg)
        @test_throws AssertionError BinaryEuropeanOption(T, Kneg)
        @test_throws AssertionError BinaryAmericanOption(T, Kneg)
        @test_throws AssertionError AsianFixedStrikeOption(T, Kneg)
        for Barrier in (BarrierOptionDownOut, BarrierOptionDownIn, BarrierOptionUpIn, BarrierOptionUpOut)
            @test_throws AssertionError Barrier(T, Kneg, D)
            @test_throws AssertionError Barrier(T, K, Kneg)
        end
        @test_throws AssertionError DoubleBarrierOption(T, K, K, Kneg)
        @test_throws AssertionError DoubleBarrierOption(T, Kneg, K, K)
        @test_throws AssertionError DoubleBarrierOption(T, K, Kneg, K)
    end

    @testset "payoffs: misc" begin
        @test FinancialMonteCarlo.maturity(Spot()) == 0.0
        @test 2.0 * EuropeanOption(T, K) == EuropeanOption(T, K) + 1.0 * EuropeanOption(T, K)
        @test 2.0 * EuropeanOption(T, K) == EuropeanOption(T, K) + EuropeanOption(T, K)
        @test_throws ErrorException payoff([1.0 1.0; 1.0 1.0], UnpricedPayoff(1.0), rfCurve, mc_anti)
    end

    @testset "smiles" begin
        K_vec = [90.0, 100.0, 110.0]
        K_vec_w = [90.0, 100.0, 110.0, 111.0]
        K_vec_w2 = [90.0, 100.0, -111.0]
        isCall_vec = [true, false, true]
        for Smile in (FinancialMonteCarlo.EuropeanOptionSmile, FinancialMonteCarlo.BinaryEuropeanOptionSmile)
            @test FinancialMonteCarlo.maturity(Smile(T, K_vec, isCall_vec)) == T
            @test_throws "K and isCall must have the same length" Smile(T, K_vec_w, isCall_vec)
            @test_throws "Strikes must be positive" Smile(T, K_vec_w2, isCall_vec)
            @test_throws "Maturity must be positive" Smile(-1.0, K_vec, isCall_vec)
        end
    end
end
