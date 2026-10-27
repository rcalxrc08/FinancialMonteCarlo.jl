include("../testutils.jl")

@testset "Statistics" begin
    model = BlackScholesProcess(sigma, Underlying(S0, d))
    eu_put = EuropeanOption(T, K, false)
    am_put = AmericanOption(T, K, false)

    @testset "variance: antithetic reduces it" begin
        @test variance(model, rfCurve, mc, eu_put) > variance(model, rfCurve, mc_anti, eu_put)
        @test all(variance(model, rfCurve, mc, [eu_put, am_put]) .>= variance(model, rfCurve, mc_anti, [eu_put, am_put]))
    end

    @testset "confidence intervals" begin
        width(ci) = ci[2] - ci[1]
        ci_std = confinter(model, rfCurve, mc, eu_put)
        ci_anti = confinter(model, rfCurve, mc_anti, eu_put, 0.99)
        @test ci_std[1] < ci_std[2]
        @test width(ci_std) > width(ci_anti)
        cis_std = confinter(model, rfCurve, mc, [eu_put, am_put])
        cis_anti = confinter(model, rfCurve, mc_anti, [eu_put, am_put], 0.99)
        @test width(cis_std[2]) > width(cis_anti[2])
    end

    @testset "terminal distribution" begin
        ST = FinancialMonteCarlo.distribution(model, rfCurve, mc, T)
        @test length(ST) == Nsim
        @test abs(exp(-r * T) * sum(ST) / length(ST) - BS_CALL_REF.fwd) < 0.8
    end
end
