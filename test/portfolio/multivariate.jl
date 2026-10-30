include("../testutils.jl")

# Pricing through a MarketDataSet/Portfolio must match pricing the model directly on the same seed.
@testset "Multivariate" begin
    bs() = BlackScholesProcess(sigma, Underlying(S0, d))
    eu_nd = EuropeanOptionND(T, K)

    @testset "bivariate Gaussian copula" begin
        model = GaussianCopulaNVariateProcess(bs(), bs(), 0.4)
        @test GaussianCopulaNVariateProcess(bs(), bs()).rho == GaussianCopulaNVariateProcess(bs(), bs(), 0.0).rho
        @test !isempty(FinancialMonteCarlo.get_parameters(model))
        @test occursin("GaussianCopulaNVariateProcess", capture_stdout(() -> display(model)))
        @test abs(pricer("ENJ_ABPL" → model, rfCurve, mc_small, "ENJ_ABPL" → eu_nd) - pricer(model, rfCurve, mc_small, eu_nd)) < 1e-8
    end

    @testset "bivariate log copula: single name legs" begin
        model = GaussianCopulaNVariateLogProcess(bs(), bs(), 0.0)
        @test GaussianCopulaNVariateLogProcess(bs(), bs()).rho == model.rho
        mkt = ("ENJ_ABPL" → model) + ("notneeded" → bs())
        portfolio = "ENJ" → 1.0 * EuropeanOption(T, K)
        portfolio += "ENJ" → 1.0 * AmericanOption(T, K)
        portfolio += "ENJ" → -1.0 * (-1.0) * BarrierOptionDownOut(T, K, D)
        @test abs(pricer(mkt, rfCurve, mc_small, portfolio) - pricer("ENJ" → bs(), rfCurve, mc_small, portfolio)) < 1e-8
    end

    @testset "trivariate $(nameof(Process))" for (Process, rho_) in [(GaussianCopulaNVariateProcess, [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]), (GaussianCopulaNVariateLogProcess, [1.0 0.0 0.1; 0.0 1.0 0.0; 0.1 0.0 1.0])]
        name = "ENJ_ABPL_TESL"
        model = Process(rho_, bs(), bs(), bs())
        @test abs(pricer(name → model, rfCurve, mc_small, name → eu_nd) - sum(pricer(model, rfCurve, mc_small, [eu_nd]))) < 1e-8
        if Process === GaussianCopulaNVariateLogProcess
            @test isfinite(pricer(name → Process(bs(), bs(), bs()), rfCurve, mc_small, name → eu_nd))
        end
    end
end
