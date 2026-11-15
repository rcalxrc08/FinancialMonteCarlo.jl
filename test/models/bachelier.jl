include("../testutils.jl")
using FinancialToolbox: normcdf, normpdf

# The Bachelier formula for S_T normal with mean the forward F and standard deviation v (FinancialToolbox has no
# Bachelier price yet: rcalxrc08/FinancialToolbox.jl#84).
function bachelier_price(F, v, K, discount, isCall)
    z = (F - K) / v
    sign_ = isCall ? 1 : -1
    return discount * (sign_ * (F - K) * normcdf(sign_ * z) + v * normpdf(z))
end

@testset "Bachelier" begin
    σ_N = 20.0
    model = BachelierProcess(σ_N, Underlying(S0, d))
    # S_T is normal, with mean the forward and standard deviation σ √T e^{(r - d) T}
    F, v = S0 * exp((r - d) * T), σ_N * sqrt(T) * exp((r - d) * T)
    ref(isCall, strike = K) = bachelier_price(F, v, strike, exp(-r * T), isCall)

    @testset "parameters" begin
        @test FinancialMonteCarlo.get_parameters(model) == [σ_N]
        @test occursin("BachelierProcess", capture_stdout(() -> display(model)))
        @test_throws AssertionError BachelierProcess(-σ_N, Underlying(S0, d))
    end

    @testset "European prices, $(nameof(typeof(config.sampling)))" for config in (mc, mc_anti, MonteCarloConfiguration(Nsim, 1; sampling = FinancialMonteCarlo.SobolMC()))
        payoffs = (fwd = Forward(T), call = EuropeanOption(T, K), put = EuropeanOption(T, K, false), otm_put = EuropeanOption(T, 80.0, false))
        refs = (fwd = S0 * exp(-d * T), call = ref(true), put = ref(false), otm_put = ref(false, 80.0))
        check_prices(model, rfCurve, config, payoffs, refs; toll = 0.6)
    end

    @testset "Sobol is accurate" begin
        mc_sobol = MonteCarloConfiguration(2^17, 1; sampling = FinancialMonteCarlo.SobolMC())
        @test abs(pricer(model, rfCurve, mc_sobol, EuropeanOption(T, K)) - ref(true)) < 1e-3
        @test abs(pricer(model, rfCurve, mc_sobol, EuropeanOption(T, K, false)) - ref(false)) < 1e-3
    end

    @testset "the discounted price with dividends is a martingale" begin
        # e^{-(r - d) t} S_t = S0 + σ W_t, and the antithetic pairs make the mean of σ W_t exactly zero
        Nstep_ = 4
        S = simulate(model, rfCurve, MonteCarloConfiguration(Nsim, Nstep_; sampling = FinancialMonteCarlo.AntitheticMC()), T)
        @test all(abs(sum(S[:, j+1]) / Nsim * exp(-(r - d) * j * T / Nstep_) - S0) < 1e-10 for j = 0:Nstep_)
        # prices can be negative
        @test any(<(0), simulate(BachelierProcess(80.0, Underlying(S0, d)), rfCurve, mc_small, T))
    end

    @testset "curves" begin
        r_curve = ZeroRate([0.0, 0.04], T)
        d_curve = FinancialMonteCarlo.Curve([0.01, 0.01], T)
        model_curve = BachelierProcess(σ_N, Underlying(S0, d_curve))
        carry = exp(FinancialMonteCarlo.integral(r_curve.r - d_curve, T))
        discount = exp(-FinancialMonteCarlo.integral(r_curve.r, T))
        mc_sobol = MonteCarloConfiguration(2^16, 3; sampling = FinancialMonteCarlo.SobolMC())
        @test abs(pricer(model_curve, r_curve, mc_sobol, EuropeanOption(T, K)) - bachelier_price(S0 * carry, σ_N * sqrt(T) * carry, K, discount, true)) < 1e-2
    end

    @testset "path dependent payoffs" begin
        check_prices(model, rfCurve, mc_small, standard_payoffs())
        check_prices(model, rfCurve, mc_small_anti, standard_put_payoffs())
    end

    @testset "Float32" begin
        model32 = BachelierProcess(20.0f0, Underlying(100.0f0, 0.01f0))
        @test pricer(model32, ZeroRate(0.02f0), mc_small, EuropeanOption(1.0f0, 100.0f0)) isa Float32
    end
end
