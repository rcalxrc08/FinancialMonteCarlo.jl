include("../testutils.jl")
using ForwardDiff
using FinancialToolbox: blsprice, normcdf, normpdf

@testset "ForwardDiff" begin
    @testset "Dual-valued volatility" begin
        model = BlackScholesProcess(ForwardDiff.Dual{Float64}(sigma, 1.0))
        # American exercise on dual paths is not checked against the reference with standard MC
        check_prices(model, und, rfCurve, mc, standard_payoffs(), Base.structdiff(BS_CALL_REF, NamedTuple{(:am,)}); toll = 0.8)
        check_prices(model, und, rfCurve, mc_anti, standard_payoffs(), BS_CALL_REF; toll = 0.6)
    end

    @testset "gradient of a European call" begin
        f(x) = pricer(BlackScholesProcess(x[1]), Underlying(x[2], d), rfCurve, mc, EuropeanOption(T, K))
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
            kou(x) = pricer(KouProcess(sigma, 5.0, x[1], x[2], x[3]), und, rfCurve, mc_jumps, EuropeanOption(T, K))
            x_kou = [0.3, 30.0, 20.0]
            @test ForwardDiff.gradient(kou, x_kou) ≈ [fd(kou, x_kou, i) for i = 1:3] rtol = 1e-4
            merton(x) = pricer(MertonProcess(sigma, 3.0, x[1], x[2]), und, rfCurve, mc_jumps, EuropeanOption(T, K))
            x_merton = [-0.02, 0.05]
            @test ForwardDiff.gradient(merton, x_merton) ≈ [fd(merton, x_merton, i) for i = 1:2] rtol = 1e-4
        end
    end

    @testset "Bachelier and shifted Black-Scholes: gradients against the formulas" begin
        # One step of a long Sobol sequence: the European price is exact up to 1e-5, and so are its derivatives
        mc_sobol = MonteCarloConfiguration(2^17, 1; sampling = FinancialMonteCarlo.SobolMC())
        # x = [S0, σ, r, d]: the Bachelier formula with mean the forward and normal volatility σ e^{(r - d) T} (not in
        # FinancialToolbox yet: rcalxrc08/FinancialToolbox.jl#84)
        function bachelier_formula(x)
            carry = exp((x[3] - x[4]) * T)
            F, v = x[1] * carry, x[2] * sqrt(T) * carry
            return exp(-x[3] * T) * ((F - K) * normcdf((F - K) / v) + v * normpdf((F - K) / v))
        end
        bachelier(x) = pricer(BachelierProcess(x[2]), Underlying(x[1], x[4]), ZeroRate(x[3]), mc_sobol, EuropeanOption(T, K))
        x_bachelier = [S0, 20.0, r, d]
        @test ForwardDiff.gradient(bachelier, x_bachelier) ≈ ForwardDiff.gradient(bachelier_formula, x_bachelier) rtol = 1e-4
        # x = [S0, σ, α]: a call on the GBM started from (1 - α) S0, with strike K - α F
        shifted_formula(x) = blsprice((1 - x[3]) * x[1], K - x[3] * x[1] * exp((r - d) * T), r, T, x[2], d)
        shifted(x) = pricer(ShiftedProcess(BlackScholesProcess(x[2]), x[3]), Underlying(x[1], d), rfCurve, mc_sobol, EuropeanOption(T, K))
        x_shifted = [S0, sigma, 0.3]
        @test ForwardDiff.gradient(shifted, x_shifted) ≈ ForwardDiff.gradient(shifted_formula, x_shifted) rtol = 1e-4
    end

    @testset "lognormal mixture: gradient against the mixture of Black-Scholes prices" begin
        # x = [S0, η₁, η₂, λ₁]; the Euler scheme leaves a small bias, in the derivatives too
        mixture_formula(x) = x[4] * blsprice(x[1], K, r, T, x[2], d) + (1 - x[4]) * blsprice(x[1], K, r, T, x[3], d)
        mixture(x) = pricer(LogNormalMixture((x[2], x[3]), (x[4],)), Underlying(x[1], d), rfCurve, MonteCarloConfiguration(2^16, 50; sampling = FinancialMonteCarlo.SobolMC()), EuropeanOption(T, K))
        x_mixture = [S0, 0.15, 0.3, 0.5]
        @test ForwardDiff.gradient(mixture, x_mixture) ≈ ForwardDiff.gradient(mixture_formula, x_mixture) rtol = 2e-2
    end

    @testset "smiles" begin
        Ks, isCall = [90.0, 100.0, 110.0], [false, true, true]
        model(s) = BlackScholesProcess(s)
        smile = FinancialMonteCarlo.EuropeanOptionSmile(T, Ks, isCall)
        vega(payoffs) = ForwardDiff.derivative(s -> pricer(model(s), und, rfCurve, mc_small, payoffs), sigma)
        @test vega(smile) ≈ vega(FinancialMonteCarlo.smile_to_options(smile))
        # the strikes are differentiated too, not only the model
        by_strike(k) = pricer(model(sigma), und, rfCurve, mc_small, FinancialMonteCarlo.EuropeanOptionSmile(T, [k, 110.0], [true, true]))
        dK = ForwardDiff.derivative(by_strike, K)
        @test dK[1] ≈ ForwardDiff.derivative(k -> pricer(model(sigma), und, rfCurve, mc_small, EuropeanOption(T, k)), K)
        @test dK[2] == 0.0
    end
end
