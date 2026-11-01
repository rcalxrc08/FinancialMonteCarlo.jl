include("../testutils.jl")
using ForwardDiff

@testset "ForwardDiff" begin
    @testset "Dual-valued volatility" begin
        model = BlackScholesProcess(ForwardDiff.Dual{Float64}(sigma, 1.0), Underlying(S0, d))
        # American exercise on dual paths is not checked against the reference with standard MC
        check_prices(model, rfCurve, mc, standard_payoffs(), Base.structdiff(BS_CALL_REF, NamedTuple{(:am,)}); toll = 0.8)
        check_prices(model, rfCurve, mc_anti, standard_payoffs(), BS_CALL_REF; toll = 0.6)
    end

    @testset "gradient of a European call" begin
        f(x) = pricer(BlackScholesProcess(x[1], Underlying(x[2], d)), rfCurve, mc, EuropeanOption(T, K))
        vega, delta_ = ForwardDiff.gradient(f, [sigma, S0])
        @test vega > 0
        @test 0 < delta_ < 1
    end

    @testset "smiles" begin
        Ks, isCall = [90.0, 100.0, 110.0], [false, true, true]
        model(s) = BlackScholesProcess(s, Underlying(S0, d))
        smile = FinancialMonteCarlo.EuropeanOptionSmile(T, Ks, isCall)
        vega(payoffs) = ForwardDiff.derivative(s -> pricer(model(s), rfCurve, mc_small, payoffs), sigma)
        @test vega(smile) ≈ vega(FinancialMonteCarlo.smile_to_options(smile))
        # the strikes are differentiated too, not only the model
        by_strike(k) = pricer(model(sigma), rfCurve, mc_small, FinancialMonteCarlo.EuropeanOptionSmile(T, [k, 110.0], [true, true]))
        dK = ForwardDiff.derivative(by_strike, K)
        @test dK[1] ≈ ForwardDiff.derivative(k -> pricer(model(sigma), rfCurve, mc_small, EuropeanOption(T, k)), K)
        @test dK[2] == 0.0
    end
end
