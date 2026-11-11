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

    @testset "jumps by inversion: pathwise derivatives in the jump law" begin
        # With antithetic and Sobol sampling the jump sizes are quantiles of uniforms, so the derivatives with respect
        # to the parameters of the jump law match finite differences on the same paths (Kou's p included). The step is
        # small, so that no path crosses the strike between the two bumps (the payoff has a kink there): with h = 1e-5
        # one path did on Julia 1.10, and dp was 4e-4 off.
        h = 1e-7
        fd(f, x, i) = (f(x .+ h .* (1:length(x) .== i)) - f(x .- h .* (1:length(x) .== i))) / (2h)
        for sampling in (FinancialMonteCarlo.AntitheticMC(), FinancialMonteCarlo.SobolMC())
            mc_jumps = MonteCarloConfiguration(2000, 5; sampling = sampling)
            kou(x) = pricer(KouProcess(sigma, 5.0, x[1], x[2], x[3], Underlying(S0, d)), rfCurve, mc_jumps, EuropeanOption(T, K))
            x_kou = [0.3, 30.0, 20.0]
            @test ForwardDiff.gradient(kou, x_kou) ≈ [fd(kou, x_kou, i) for i = 1:3] rtol = 1e-4
            merton(x) = pricer(MertonProcess(sigma, 3.0, x[1], x[2], Underlying(S0, d)), rfCurve, mc_jumps, EuropeanOption(T, K))
            x_merton = [-0.02, 0.05]
            @test ForwardDiff.gradient(merton, x_merton) ≈ [fd(merton, x_merton, i) for i = 1:2] rtol = 1e-4
        end
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
