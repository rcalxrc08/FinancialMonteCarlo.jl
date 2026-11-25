include("../testutils.jl")
using Statistics: var, std
using Distributions: TDist, quantile

@testset "Statistics" begin
    model = BlackScholesProcess(sigma)
    eu_put = EuropeanOption(T, K, false)
    am_put = AmericanOption(T, K, false)

    @testset "variance: antithetic reduces it" begin
        @test variance(model, und, rfCurve, mc, eu_put) > variance(model, und, rfCurve, mc_anti, eu_put)
        @test all(variance(model, und, rfCurve, mc, [eu_put, am_put]) .>= variance(model, und, rfCurve, mc_anti, [eu_put, am_put]))
    end

    @testset "confidence intervals" begin
        width(ci) = ci[2] - ci[1]
        ci_std = confinter(model, und, rfCurve, mc, eu_put)
        ci_anti = confinter(model, und, rfCurve, mc_anti, eu_put, 0.99)
        @test ci_std[1] < ci_std[2]
        @test width(ci_std) > width(ci_anti)
        cis_std = confinter(model, und, rfCurve, mc, [eu_put, am_put])
        cis_anti = confinter(model, und, rfCurve, mc_anti, [eu_put, am_put], 0.99)
        @test width(cis_std[2]) > width(cis_anti[2])
    end

    @testset "antithetic pairs are the independent samples (#76)" begin
        # the variance per path is twice the sample variance of the averages of the pairs (path i, path i + Nsim / 2)
        cfg = MonteCarloConfiguration(2000, 1; sampling = FinancialMonteCarlo.AntitheticMC(), seed = 3)
        straddle = EuropeanOption(T, K) + EuropeanOption(T, K, false)
        FinancialMonteCarlo.set_seed!(cfg)
        values = payoff(simulate(model, und, rfCurve, cfg, T), straddle, rfCurve, cfg, T)
        pairs = (values[1:1000] .+ values[1001:2000]) ./ 2
        @test variance(model, und, rfCurve, cfg, straddle) ≈ 2 * var(pairs)
        @test variance(model, und, rfCurve, cfg, [eu_put, EuropeanOption(T, K)]) ≈ [variance(model, und, rfCurve, cfg, eu_put), variance(model, und, rfCurve, cfg, EuropeanOption(T, K))]
        # the interval of the mean of 1000 pairs: 999 degrees of freedom
        lo, hi = confinter(model, und, rfCurve, cfg, straddle, 0.95)
        @test (hi - lo) / 2 ≈ quantile(TDist(999), 0.975) * std(pairs) / sqrt(1000)
        # an even payoff: the pairs are positively correlated, antithetic sampling increases the variance (it was
        # reported 30% too small); a monotone one: it decreases it
        mc_std, mc_pairs = MonteCarloConfiguration(20_000, 1; seed = 1), MonteCarloConfiguration(20_000, 1; sampling = FinancialMonteCarlo.AntitheticMC(), seed = 1)
        @test variance(model, und, rfCurve, mc_pairs, straddle) > 1.5 * variance(model, und, rfCurve, mc_std, straddle)
        @test variance(model, und, rfCurve, mc_pairs, EuropeanOption(T, K)) < variance(model, und, rfCurve, mc_std, EuropeanOption(T, K))
        # the batches pool the moments of their pairs
        batched = MonteCarloConfiguration(20_000, 1; sampling = FinancialMonteCarlo.AntitheticMC(), mode = FinancialMonteCarlo.MultiThreading(4), seed = 1)
        @test variance(model, und, rfCurve, batched, straddle) ≈ variance(model, und, rfCurve, mc_pairs, straddle) rtol = 0.1
    end

    # A vector of payoffs of different types has an abstract element type: the result must still be inferred
    @testset "inferred for a vector of payoffs" begin
        payoffs = [eu_put, am_put]
        @inferred pricer(model, und, rfCurve, mc_small, payoffs)
        @inferred variance(model, und, rfCurve, mc_small, payoffs)
        @inferred confinter(model, und, rfCurve, mc_small, payoffs)
        @inferred delta(model, und, rfCurve, mc_small, payoffs)
        @inferred rho(model, und, rfCurve, mc_small, payoffs)
    end

    @testset "the metrics stay in the precision of the inputs" begin
        model_32 = BlackScholesProcess(Float32(sigma))
        curve_32 = ZeroRate(Float32(r))
        eu_32 = EuropeanOption(Float32(T), Float32(K))
        smile_32 = FinancialMonteCarlo.EuropeanOptionSmile(Float32(T), Float32[K, K+10], [true, false])
        threaded = MonteCarloConfiguration(Nsim_small, Nstep_small; mode = FinancialMonteCarlo.MultiThreading(2), seed = [1, 2])
        for config in (mc_small, threaded)
            @test pricer(model_32, Underlying(Float32(S0), Float32(d)), curve_32, config, eu_32) isa Float32
            @test variance(model_32, Underlying(Float32(S0), Float32(d)), curve_32, config, eu_32) isa Float32
            @test confinter(model_32, Underlying(Float32(S0), Float32(d)), curve_32, config, eu_32) isa Tuple{Float32, Float32}
            @test confinter(model_32, Underlying(Float32(S0), Float32(d)), curve_32, config, [eu_32, eu_32]) isa Vector{Tuple{Float32, Float32}}
            @test confinter(model_32, Underlying(Float32(S0), Float32(d)), curve_32, config, smile_32) isa Vector{Tuple{Float32, Float32}}
        end
        # no Float64 constant of the package promotes them: default dividend, jump drift, basket payoff, position weights
        und_32 = Underlying(Float32(S0), Float32(d))
        @test pricer(BlackScholesProcess(Float32(sigma)), Underlying(Float32(S0)), curve_32, mc_small, eu_32) isa Float32
        @test pricer(MertonProcess(Float32(sigma), 3.0f0, 0.03f0, 0.02f0), und_32, curve_32, mc_small, eu_32) isa Float32
        @test pricer(MertonProcess(Float32(sigma), 3.0f0, 0.03f0, 0.02f0), und_32, curve_32, mc_small_anti, eu_32) isa Float32
        @test pricer(KouProcess(Float32(sigma), 3.0f0, 0.3f0, 30.0f0, 20.0f0), und_32, curve_32, MonteCarloConfiguration(Nsim_small, Nstep_small; sampling = FinancialMonteCarlo.SobolMC()), eu_32) isa Float32
        @test pricer(GaussianCopulaNVariateProcess(model_32, model_32, 0.4f0), (und_32, und_32), curve_32, mc_small, EuropeanOptionND(Float32(T), Float32(K))) isa Float32
        @test pricer(model_32, Underlying(Float32(S0), Float32(d)), curve_32, mc_small, eu_32 + Forward(Float32(T))) isa Float32
        @test pricer("A" → (model_32, und_32), curve_32, mc_small, "A" → eu_32) isa Float32
        @test pricer(model_32, Underlying(Float32(S0), Float32(d)), curve_32, mc_small, eu_32 + Spot()) isa Float32 # #39
        @test FinancialMonteCarlo.integral(ZeroRate(FinancialMonteCarlo.Curve(Float32[r, 2r], Float32(T))).r, 0.0f0) isa Float32
    end

    @testset "a payoff of a vector is the payoff alone" begin
        @test variance(model, und, rfCurve, mc_small, [eu_put, am_put])[1] ≈ variance(model, und, rfCurve, mc_small, eu_put)
        @test all(confinter(model, und, rfCurve, mc_small, [eu_put, am_put])[1] .≈ confinter(model, und, rfCurve, mc_small, eu_put))
    end

    @testset "terminal distribution" begin
        ST = FinancialMonteCarlo.distribution(model, und, rfCurve, mc, T)
        @test length(ST) == Nsim
        @test abs(exp(-r * T) * sum(ST) / length(ST) - BS_CALL_REF.fwd) < 0.8
    end
end
