include("../testutils.jl")
using FinancialToolbox: blsprice

# The lognormal mixture of Brigo and Mercurio: the law of S_t is the mixture of the lognormal laws, so European prices are
# the weighted sums of the Black-Scholes prices (the closed form of FinancialCalibrator's AnalyticMethod).
@testset "LogNormal mixture" begin
    und = Underlying(S0, d)
    η, λ = (0.15, 0.3), (0.5,)
    model = LogNormalMixture(η, λ)
    mixture_price(η, λ, K, isCall, T_ = T; spot = S0, strike_shift = 0.0) = sum(w * blsprice(spot, K - strike_shift, r, T_, η_i, d, isCall) for (w, η_i) in zip((λ..., 1 - sum(λ)), η))
    mc_sobol = MonteCarloConfiguration(2^16, 50; sampling = FinancialMonteCarlo.SobolMC())

    @testset "parameters" begin
        @test FinancialMonteCarlo.get_parameters(model) == [0.15, 0.3, 0.5]
        @test FinancialMonteCarlo.set_parameters!(model, [0.2, 0.25, 0.3]).η == (0.2, 0.25)
        @test occursin("LogNormalMixture", capture_stdout(() -> display(model)))
        @test_throws DomainError LogNormalMixture((-0.2, 0.2, 0.2), (0.1, 0.1))
        @test_throws DomainError LogNormalMixture((0.2, 0.2, 0.2), (-0.1, 0.1))
        @test_throws DomainError LogNormalMixture((0.2, 0.2, 0.2), (1.0, 1.0))
        @test_throws DimensionMismatch LogNormalMixture((0.2, 0.2, 0.2), (0.2, 0.2, 0.2))
    end

    @testset "European prices against the mixture of Black-Scholes prices" begin
        for (η_, λ_) in ((η, λ), ((0.2, 0.3, 0.25), (0.3, 0.4))), K_ in (70.0, 100.0, 130.0), isCall in (true, false)
            m = LogNormalMixture(η_, λ_)
            @test abs(pricer(m, und, rfCurve, mc_sobol, EuropeanOption(T, K_, isCall)) - mixture_price(η_, λ_, K_, isCall)) < 0.05
        end
        # the law is the mixture at every time, not only at the end of the paths
        @test abs(pricer(model, und, rfCurve, mc_sobol, EuropeanOption(T / 4, K)) - mixture_price(η, λ, K, true, T / 4)) < 0.05
        @test abs(pricer(model, und, rfCurve, mc_sobol, Forward(T)) - S0 * exp(-d * T)) < 0.02
    end

    @testset "$(nameof(typeof(config.sampling)))" for config in (mc, mc_anti)
        refs = (fwd = S0 * exp(-d * T), eu = mixture_price(η, λ, K, true))
        check_prices(model, und, rfCurve, config, standard_payoffs(), refs; toll = 0.6)
    end

    @testset "equal volatilities: Black-Scholes" begin
        flat = LogNormalMixture((sigma, sigma), (0.3,))
        check_prices(flat, und, rfCurve, mc, standard_payoffs(), BS_CALL_REF; toll = 0.8)
        @test abs(pricer(flat, und, rfCurve, mc_sobol, EuropeanOption(T, K)) - BS_CALL_REF.eu) < 0.02
    end

    @testset "shifted, against the shifted closed form" begin
        α = 0.2
        shifted = ShiftedProcess(model, α)
        F = S0 * exp((r - d) * T)
        for K_ in (80.0, 100.0, 130.0)
            @test abs(pricer(shifted, und, rfCurve, mc_sobol, EuropeanOption(T, K_)) - mixture_price(η, λ, K_, true; spot = (1 - α) * S0, strike_shift = α * F)) < 0.05
        end
    end

    @testset "curves" begin
        r_curve = ZeroRate(FinancialMonteCarlo.Curve([0.0, 0.04], T))
        model_curve = LogNormalMixture(η, λ)
        @test abs(pricer(model_curve, Underlying(S0, FinancialMonteCarlo.Curve([0.01, 0.01], T)), r_curve, mc_sobol, Forward(T)) - S0 * exp(-d * T)) < 0.02
        check_prices(model_curve, Underlying(S0, FinancialMonteCarlo.Curve([0.01, 0.01], T)), r_curve, mc_small_anti, standard_payoffs())
    end

    @testset "Float32" begin
        model32 = LogNormalMixture((0.15f0, 0.3f0), (0.5f0,))
        @test pricer(model32, Underlying(100.0f0, 0.01f0), ZeroRate(0.02f0), mc_small, EuropeanOption(1.0f0, 100.0f0)) isa Float32
    end
end
