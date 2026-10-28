include("../testutils.jl")

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

    @testset "parameter validation" begin
        @test_throws AssertionError KouProcess(-sigma, lam, p, lamp, lamm, Underlying(S0, d))
        @test_throws AssertionError KouProcess(sigma, -lam, p, lamp, lamm, Underlying(S0, d))
        @test_throws AssertionError KouProcess(sigma, lam, -p, lamp, lamm, Underlying(S0, d))
        @test_throws AssertionError KouProcess(sigma, lam, p, -lamp, lamm, Underlying(S0, d))
        @test_throws AssertionError KouProcess(sigma, lam, p, 0.5, lamm, Underlying(S0, d))
        @test_throws AssertionError KouProcess(sigma, lam, p, lamp, -lamm, Underlying(S0, d))
    end
end
