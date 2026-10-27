include("../testutils.jl")

const MERTON_REF = (fwd = 99.1188767166039, eu = 9.084327245917533, barrier = 7.880881290426765, asian_float = 5.129020349580892)

@testset "Merton" begin
    lam, mu1, sigma1 = 5.0, 0.03, 0.02
    model = MertonProcess(sigma, lam, mu1, sigma1, Underlying(S0, d))

    @testset "standard" begin
        check_prices(model, rfCurve, mc, standard_payoffs(), MERTON_REF; toll = 0.8)
    end

    @testset "antithetic" begin
        check_prices(model, rfCurve, mc_anti, standard_payoffs(), MERTON_REF; toll = 0.6)
    end

    @testset "parameter validation" begin
        @test_throws AssertionError MertonProcess(-sigma, lam, mu1, sigma1, Underlying(S0, d))
        @test_throws AssertionError MertonProcess(sigma, lam, mu1, -sigma1, Underlying(S0, d))
        @test_throws AssertionError MertonProcess(sigma, -lam, mu1, sigma1, Underlying(S0, d))
    end
end
