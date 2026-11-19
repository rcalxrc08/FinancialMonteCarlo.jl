include("testutils.jl")

# A payoff with no `payoff` method, to check the generic fallback errors out.
# An adaptor that negates the strikes, to check Adapt.adapt_structure of a smile does not validate it again (#73).
struct NegatedStrikes end
FinancialMonteCarlo.Adapt.adapt_storage(::NegatedStrikes, x::AbstractVector{Float64}) = -x

struct UnpricedPayoff <: FinancialMonteCarlo.AbstractPayoff{Float64}
    T::Float64
end

@testset "Structs" begin
    Kneg = -K
    Tneg = -T

    @testset "MonteCarloConfiguration" begin
        FMC = FinancialMonteCarlo
        @test MonteCarloConfiguration(Nsim, Nstep; seed = 1).rng.seed == 1
        @test MonteCarloConfiguration(Nsim, Nstep).mode isa FMC.SerialMode
        @test_throws DomainError MonteCarloConfiguration(-Nsim, Nstep)
        @test_throws DomainError MonteCarloConfiguration(Nsim, -Nstep)
        @test_throws ArgumentError MonteCarloConfiguration(Nsim + 1, Nstep; sampling = FMC.AntitheticMC())
        @test_throws ArgumentError MonteCarloConfiguration(Nsim, Nstep; sampling = FMC.SobolMC(), mode = FMC.CudaMode())
        # the mode is part of the method
        mc_mt = MonteCarloConfiguration(Nsim, Nstep; seed = 3, mode = FMC.MultiThreading(2))
        @test FMC.execution_mode(mc_mt) === FMC.MultiThreading(2)
        mc_serial = FMC.with_mode(mc_mt, FMC.SerialMode())
        @test FMC.execution_mode(mc_serial) === FMC.SerialMode()
        @test (mc_serial.Nsim, mc_serial.Nstep, mc_serial.rng.seed) == (Nsim, Nstep, 3)
        @test [state.seed for state in FMC.with_mode(mc_serial, FMC.MultiThreading(2)).rng] == [4, 5]
        # the names and positional constructors of version 0.5 are gone (#77)
        @test !isdefined(FMC, :SobolMode) && !isdefined(FMC, :BaseMode) && !isdefined(FMC, :AbstractMonteCarloMethod)
        @test_throws MethodError MonteCarloConfiguration(Nsim, Nstep, FMC.AntitheticMC())
    end

    @testset "Underlying and curves" begin
        rf_curve = ZeroRate([0.00, 0.02], T)
        @test FinancialMonteCarlo.Curve([0.0, 0.1], [0.0, 0.1]) isa FinancialMonteCarlo.Curve{2, Float64, Float64}
        @test_throws DomainError Underlying(-S0)
        @test_throws DomainError Underlying(-S0, rf_curve.r)

        # trapezoidal inside the nodes, flat outside them
        c = ZeroRate([0.01, 0.02, 0.05], 1.0).r # nodes at 0, 0.5, 1
        @test FinancialMonteCarlo.integral(c, 1.0) ≈ 0.025
        @test FinancialMonteCarlo.integral(c, 2.0) ≈ 0.025 + 0.05
        c_late = FinancialMonteCarlo.Curve([0.03, 0.05], [0.5, 1.0])
        @test FinancialMonteCarlo.integral(c_late, 0.25) ≈ 0.03 * 0.25
        @test FinancialMonteCarlo.integral(c_late, 1.0) ≈ 0.03 * 0.5 + 0.04 * 0.5
        @test FinancialMonteCarlo.integral(FinancialMonteCarlo.Curve([0.03], [1.0]), 2.0) ≈ 0.06

        # arithmetic with numbers on the values; between curves only + and -, lazy and exact (#81)
        Curve = FinancialMonteCarlo.Curve
        c1, c2, c3 = Curve([0.01, 0.02], 1.0), Curve([0.03, 0.03, 0.03], 1.0), Curve([0.01, 0.02], [0.5, 1.0])
        @test (c1 + 1).values ≈ [1.01, 1.02] && (1 - c1).values ≈ [0.99, 0.98] && (2 * c1).values ≈ [0.02, 0.04] && (c1 / 2).values ≈ [0.005, 0.01]
        @test c3 - c2 isa FinancialMonteCarlo.CurveDifference && c3 + c2 isa FinancialMonteCarlo.CurveDifference
        for t in (0.0, 0.25, 0.5, 0.7, 1.0, 1.5)
            @test FinancialMonteCarlo.integral(c3 - c2, t) ≈ FinancialMonteCarlo.integral(c3, t) - FinancialMonteCarlo.integral(c2, t)
            @test FinancialMonteCarlo.integral(c3 + c2, t) ≈ FinancialMonteCarlo.integral(c3, t) + FinancialMonteCarlo.integral(c2, t)
            @test FinancialMonteCarlo.integral(c3 - c2 - 0.01, t) ≈ FinancialMonteCarlo.integral(c3, t) - FinancialMonteCarlo.integral(c2, t) - 0.01 * t
            @test FinancialMonteCarlo.integral(0.01 + (c3 - c2), t) ≈ FinancialMonteCarlo.integral(c3, t) - FinancialMonteCarlo.integral(c2, t) + 0.01 * t
            @test FinancialMonteCarlo.integral(-(c3 - c2), t) ≈ FinancialMonteCarlo.integral(c2, t) - FinancialMonteCarlo.integral(c3, t)
            @test FinancialMonteCarlo.value_at(c1 - c3, t) ≈ FinancialMonteCarlo.value_at(c1, t) - FinancialMonteCarlo.value_at(c3, t)
        end
        @test_throws MethodError c1 * c2
        @test_throws MethodError c1 / c2
        # a drift r - d with a dividend curve finer than the rate curve keeps its sign
        model = BlackScholesProcess(sigma, Underlying(S0, Curve([0.0, 0.0, 0.0], T)))
        @test pricer(model, ZeroRate([0.05, 0.05], T), mc_small, Forward(T)) ≈ S0 rtol = 0.02
        @test_throws ArgumentError Curve([0.01, 0.02], [1.0, 0.5])
        # a type of the package, not a Dictionary with arithmetic defined on it (#42)
        @test parentmodule(Curve) === FinancialMonteCarlo
        @test FinancialMonteCarlo.ImpliedCurve([0.01, 0.02], [1.0, 2.0]) isa Curve{3, Float64, Float64}
        # isbits (#79), and inferred from tuples
        @test isbits(BlackScholesProcess(sigma, Underlying(S0, c3))) && isbits(ZeroRate([0.01, 0.02], T))
        @test (@inferred FinancialMonteCarlo.ImpliedCurve((0.01, 0.02), (1.0, 2.0))) == FinancialMonteCarlo.ImpliedCurve([0.01, 0.02], [1.0, 2.0])
        @test (@inferred Curve((0.01, 0.02), (0.5, 1.0))) == c3
        @test (@inferred ZeroRate((0.01, 0.02), 1.0)).r == c1
        @test FinancialMonteCarlo.integral(Curve((0.01,), 1.0), 2.0) ≈ 0.02
        @test FinancialMonteCarlo.integral(FinancialMonteCarlo.ImpliedCurve([0.01, 0.02], [1.0, 2.0]), 2.0) ≈ 0.04
    end

    @testset "market data and portfolios" begin
        model = GaussianCopulaNVariateProcess(BlackScholesProcess(sigma, Underlying(S0, d)), BlackScholesProcess(sigma, Underlying(S0, d)), 0.0)
        mkt = ("eni", "aap") → model
        @test_throws ArgumentError mkt + mkt
        @test_throws ArgumentError mkt + ("aap" → BlackScholesProcess(sigma, Underlying(S0, d))) # a name appears once
        @test_throws ArgumentError "ENJ" → model # a multivariate model needs one name per component
        @test_throws ArgumentError ("a", "b", "c") → model
        @test_throws ArgumentError ("a", "b") → Forward(1.0) # a single name payoff needs a single name
        @test_throws ArgumentError "ENI" → EuropeanOptionND(T, K) # a basket payoff needs a tuple of names (#36)
        # names joined by `_` are deprecated: they become the tuple of the names
        @test (@test_deprecated "eni_aap" → model) == mkt
        @test_throws ArgumentError "c_i" → Forward(1.0) # deprecated, and the tuple of two names is not a single name
    end

    @testset "payoffs: negative maturity" begin
        @test_throws DomainError Forward(Tneg)
        @test_throws DomainError EuropeanOption(Tneg, K)
        @test_throws DomainError EuropeanOptionND(Tneg, K)
        @test_throws DomainError AmericanOption(Tneg, K)
        @test_throws DomainError BermudanOption([Tneg], K)
        @test_throws DomainError BinaryAmericanOption(Tneg, K)
        @test_throws DomainError BinaryEuropeanOption(Tneg, K)
        @test_throws DomainError BarrierOptionDownOut(Tneg, K, D)
        @test_throws DomainError BarrierOptionUpOut(Tneg, K, D)
        @test_throws DomainError BarrierOptionUpIn(Tneg, K, D)
        @test_throws DomainError BarrierOptionDownIn(Tneg, K, D)
        @test_throws DomainError AsianFloatingStrikeOption(Tneg)
        @test_throws DomainError AsianFixedStrikeOption(Tneg, K)
        @test_throws DomainError DoubleBarrierOption(Tneg, K, K * 10, D)
    end

    @testset "payoffs: negative strike or barrier" begin
        @test_throws DomainError EuropeanOption(T, Kneg)
        @test_throws DomainError EuropeanOptionND(T, Kneg)
        @test_throws DomainError AmericanOption(T, Kneg)
        @test_throws DomainError BermudanOption([T], Kneg)
        @test_throws DomainError BinaryEuropeanOption(T, Kneg)
        @test_throws DomainError BinaryAmericanOption(T, Kneg)
        @test_throws DomainError AsianFixedStrikeOption(T, Kneg)
        for Barrier in (BarrierOptionDownOut, BarrierOptionDownIn, BarrierOptionUpIn, BarrierOptionUpOut)
            @test_throws DomainError Barrier(T, Kneg, D)
            @test_throws DomainError Barrier(T, K, Kneg)
        end
        @test_throws DomainError DoubleBarrierOption(T, K, K, Kneg)
        @test_throws DomainError DoubleBarrierOption(T, Kneg, K, K)
        @test_throws DomainError DoubleBarrierOption(T, K, Kneg, K)
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
            # strikes and flags are vectors (#37): a range or a view works, a matrix does not
            smile_range = Smile(T, 90.0:10.0:110.0, view(isCall_vec, :))
            S = simulate(BlackScholesProcess(sigma, Underlying(S0, d)), rfCurve, mc_small, T)
            @test payoff(S, smile_range, rfCurve, mc_small) == payoff(S, Smile(T, K_vec, isCall_vec), rfCurve, mc_small)
            @test_throws MethodError Smile(T, reshape(K_vec, 1, :), reshape(isCall_vec, 1, :))
            # to_mode moves the arrays without checking the smile again (#73): an adaptor that negates the strikes passes
            smile = Smile(T, K_vec, isCall_vec)
            @test FinancialMonteCarlo.to_mode(smile, FinancialMonteCarlo.SerialMode()).K == K_vec
            @test FinancialMonteCarlo.Adapt.adapt(NegatedStrikes(), smile).K == -K_vec
        end
    end

    @testset "one-dimensional inputs are vectors (#37)" begin
        @test BermudanOption(0.5:0.25:1.0, K).T_ex == [0.5, 0.75, 1.0]
        @test_throws MethodError BermudanOption([0.5 1.0], K)
        @test FinancialMonteCarlo.integral(ZeroRate(0.01:0.01:0.02, T).r, T) ≈ FinancialMonteCarlo.integral(ZeroRate([0.01, 0.02], T).r, T)
        @test_throws MethodError ZeroRate([0.01 0.02], T)
        model = BlackScholesProcess(sigma, Underlying(S0, d))
        @test FinancialMonteCarlo.set_parameters!(model, view([0.3, 1.0], 1:1)).σ == 0.3
    end
end
