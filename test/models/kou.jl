include("../testutils.jl")
using Distributions

const KOU_REF = (fwd = 99.41332633109904, eu = 10.347332240535199, barrier = 8.860123655599818, asian_float = 5.81798437145069)

@testset "Kou" begin
    p, lam, lamp, lamm = 0.3, 5.0, 30.0, 20.0
    model = KouProcess(sigma, lam, p, lamp, lamm, Underlying(S0, d))

    @testset "standard" begin
        check_prices(model, rfCurve, mc, standard_payoffs(), KOU_REF; toll = 0.8)
    end

    @testset "antithetic" begin
        check_prices(model, rfCurve, mc_anti, standard_payoffs(), KOU_REF; toll = 0.8)
    end

    @testset "sobol" begin
        mc_sobol = MonteCarloConfiguration(Nsim, Nstep; sampling = FinancialMonteCarlo.SobolMC())
        check_prices(model, rfCurve, mc_sobol, standard_payoffs(), KOU_REF; toll = 0.8)
        # the jumps are quasi random too: the forward is close to its exact value
        @test abs(pricer(model, rfCurve, mc_sobol, Forward(T)) - S0 * exp(-d * T)) < 0.05
    end

    @testset "jump sizes by inversion" begin
        # double exponential: downwards with probability 1 - p, continuous at u = 1 - p
        q = 1 - p
        for u in (0.01, 0.2, 0.5, 0.9, 0.999)
            expected = u < q ? quantile(Exponential(1 / lamm), (q - u) / q) : quantile(Exponential(1 / lamp), (u - q) / p)
            @test FinancialMonteCarlo.jump_quantile(model, u) ≈ (u < q ? -expected : expected)
        end
        @test abs(FinancialMonteCarlo.jump_quantile(model, prevfloat(q))) < 1e-12
        @test abs(FinancialMonteCarlo.jump_quantile(model, q)) < 1e-12
    end

    @testset "parameter validation" begin
        @test_throws AssertionError KouProcess(-sigma, lam, p, lamp, lamm, Underlying(S0, d))
        @test_throws AssertionError KouProcess(sigma, -lam, p, lamp, lamm, Underlying(S0, d))
        @test_throws AssertionError KouProcess(sigma, lam, -p, lamp, lamm, Underlying(S0, d))
        @test_throws AssertionError KouProcess(sigma, lam, p, -lamp, lamm, Underlying(S0, d))
        @test_throws AssertionError KouProcess(sigma, lam, p, 0.5, lamm, Underlying(S0, d))
        @test_throws AssertionError KouProcess(sigma, lam, p, lamp, -lamm, Underlying(S0, d))
    end
end
