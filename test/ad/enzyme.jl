include("../testutils.jl")
using Enzyme, ForwardDiff, Random

@testset "Enzyme" begin
    @testset "Kou European call: gradient matches ForwardDiff" begin
        # The default MersenneTwister: its scalar draws (jump times and sizes) are inactive in the extension.
        # x = [S0, σ, λ, p, λ₊, λ₋, r, d, T]
        x0 = [S0, sigma, 5.0, 0.3, 30.0, 20.0, r, d, T]
        mc_kou = MonteCarloConfiguration(100, Nstep_small)
        f(x) = pricer(KouProcess(x[2], x[3], x[4], x[5], x[6], Underlying(x[1], x[8])), ZeroRate(x[7]), mc_kou, EuropeanOption(x[9], K))
        g_ref = ForwardDiff.gradient(f, x0)
        g = Enzyme.gradient(Enzyme.Reverse, Enzyme.Const(f), x0)[1]
        @test g ≈ g_ref rtol = 1e-8
        @test collect(Enzyme.gradient(Enzyme.Forward, Enzyme.Const(f), x0)[1]) ≈ g_ref rtol = 1e-8
        @test 0 < g[1] < 1 # delta
        @test g[2] > 0 # vega
    end

    @testset "two-leg position built inside the function (#36, #40)" begin
        mc_xo = MonteCarloConfiguration(200, Nstep_small; rng = Xoshiro())
        f(x) = pricer(BlackScholesProcess(x[2], Underlying(x[1], d)), ZeroRate(r), mc_xo, 2 * EuropeanOption(T, x[3]) - EuropeanOption(T, x[3], false))
        x0 = [S0, sigma, K]
        # on Julia 1.10 Enzyme cannot prove statically that the integer weights stored next to the active payoffs are
        # inactive, and asks for runtime activity
        mode = VERSION < v"1.11" ? Enzyme.set_runtime_activity(Enzyme.Reverse) : Enzyme.Reverse
        @test Enzyme.gradient(mode, Enzyme.Const(f), x0)[1] ≈ ForwardDiff.gradient(f, x0) rtol = 1e-10
    end

    @testset "the rules of the kernels: reverse and forward mode match ForwardDiff" begin
        # every kernel with a rule in the extension (bm_step!, exp_paths!, affine_carry!, add_carry!, heston_step!,
        # lnm_step!, add_jump_sizes!, path_range), antithetic, with a rate curve and three payoffs on the same paths
        x0 = [S0, sigma, r, d]
        mc_anti_xo = MonteCarloConfiguration(200, Nstep_small; sampling = FinancialMonteCarlo.AntitheticMC(), rng = Xoshiro())
        # the three payoffs priced one by one: a position of payoffs of different types holds them in an abstract vector,
        # which Enzyme does not differentiate (#65)
        legs = (EuropeanOption(T, K), AsianFixedStrikeOption(T, K), BarrierOptionDownOut(T, K, 85.0))
        models = (x -> BlackScholesProcess(x[2], Underlying(x[1], x[4])), x -> BachelierProcess(100 * x[2], Underlying(x[1], x[4])), x -> HestonProcess(0.5, x[2], 0.0, 2.0, -0.7, 0.04, Underlying(x[1], x[4])), x -> KouProcess(x[2], 3.0, 0.3, 30.0, 20.0, Underlying(x[1], x[4])), x -> MertonProcess(x[2], 3.0, -0.05, 0.1, Underlying(x[1], x[4])), x -> LogNormalMixture((x[2], 0.3), (0.5,), Underlying(x[1], x[4])), x -> ShiftedProcess(BlackScholesProcess(x[2], Underlying(x[1], x[4])), 0.2))
        @testset "model $k" for (k, build) in enumerate(models)
            f(x) = sum(leg -> pricer(build(x), ZeroRate([x[3], 2 * x[3]], T), mc_anti_xo, leg), legs)
            g_ref = ForwardDiff.gradient(f, x0)
            @test Enzyme.gradient(Enzyme.Reverse, Enzyme.Const(f), x0)[1] ≈ g_ref rtol = 1e-10
            @test collect(Enzyme.gradient(Enzyme.Forward, Enzyme.Const(f), x0)[1]) ≈ g_ref rtol = 1e-10
        end
    end

    @testset "models on a dividend curve, without runtime activity (#79)" begin
        # the nodes of a Curve are SVectors: an UnderlyingVec is isbits, so a model built from active parameters on it
        # holds no constant heap memory next to the active numbers (EnzymeRuntimeActivityError with Vector nodes)
        underlying = Underlying(S0, FinancialMonteCarlo.ImpliedCurve([0.01, 0.012, 0.015], [0.5, 1.0, 2.0]))
        @noinline escape(m) = m.σ * m.underlying.S0 # the model escapes into a call, as in every pricer
        f_one(x) = escape(BlackScholesProcess(x[1], underlying))
        @test Enzyme.gradient(Enzyme.Reverse, Enzyme.Const(f_one), [sigma])[1] ≈ [S0]
        # nested models rebuilt by set_parameters!, as FinancialCalibrator's vol_errors does
        shifted = ShiftedProcess(LogNormalMixture((0.2, 0.3), (0.6,), underlying), 0.1)
        @noinline escape_mixture(m) = m.model.η[1] * m.underlying.S0 + m.model.λ[1] * m.α
        f_nested(x) = escape_mixture(FinancialMonteCarlo.set_parameters!(shifted, x))
        x0 = FinancialMonteCarlo.get_parameters(shifted)
        @test Enzyme.gradient(Enzyme.Reverse, Enzyme.Const(f_nested), x0)[1] ≈ ForwardDiff.gradient(f_nested, x0) rtol = 1e-12
        # a price on the curve
        mc_xo = MonteCarloConfiguration(500, Nstep_small; rng = Xoshiro())
        f_price(x) = pricer(BlackScholesProcess(x[2], Underlying(x[1], underlying.d)), ZeroRate(x[3]), mc_xo, EuropeanOption(T, K))
        x1 = [S0, sigma, r]
        @test Enzyme.gradient(Enzyme.Reverse, Enzyme.Const(f_price), x1)[1] ≈ ForwardDiff.gradient(f_price, x1) rtol = 1e-10
    end

    @testset "Black-Scholes, Xoshiro" begin
        x0 = [S0, sigma, r, d]
        mc_xo = MonteCarloConfiguration(2000, Nstep_small; rng = Xoshiro())
        f(x) = pricer(BlackScholesProcess(x[2], Underlying(x[1], x[4])), ZeroRate(x[3]), mc_xo, EuropeanOption(T, K))
        g = Enzyme.gradient(Enzyme.Reverse, Enzyme.Const(f), x0)[1]
        @test g ≈ ForwardDiff.gradient(f, x0) rtol = 1e-10
    end
end
