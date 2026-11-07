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

    @testset "the metrics stay in the precision of the inputs" begin
        model_32 = BlackScholesProcess(Float32(sigma), Underlying(Float32(S0), Float32(d)))
        curve_32 = ZeroRate(Float32(r))
        eu_32 = EuropeanOption(Float32(T), Float32(K))
        smile_32 = FinancialMonteCarlo.EuropeanOptionSmile(Float32(T), Float32[K, K + 10], [true, false])
        threaded = MonteCarloConfiguration(Nsim_small, Nstep_small; mode = FinancialMonteCarlo.MultiThreading(2), seed = [1, 2])
        for config in (mc_small, threaded)
            @test pricer(model_32, curve_32, config, eu_32) isa Float32
            @test variance(model_32, curve_32, config, eu_32) isa Float32
            @test confinter(model_32, curve_32, config, eu_32) isa Tuple{Float32, Float32}
            @test confinter(model_32, curve_32, config, [eu_32, eu_32]) isa Vector{Tuple{Float32, Float32}}
            @test confinter(model_32, curve_32, config, smile_32) isa Vector{Tuple{Float32, Float32}}
        end
    end

    @testset "a payoff of a vector is the payoff alone" begin
        @test variance(model, rfCurve, mc_small, [eu_put, am_put])[1] ≈ variance(model, rfCurve, mc_small, eu_put)
        @test all(confinter(model, rfCurve, mc_small, [eu_put, am_put])[1] .≈ confinter(model, rfCurve, mc_small, eu_put))
    end

    @testset "terminal distribution" begin
        ST = FinancialMonteCarlo.distribution(model, rfCurve, mc, T)
        @test length(ST) == Nsim
        @test abs(exp(-r * T) * sum(ST) / length(ST) - BS_CALL_REF.fwd) < 0.8
    end
end
