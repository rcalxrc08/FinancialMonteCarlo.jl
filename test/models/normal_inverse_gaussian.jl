include("../testutils.jl")

const NIG_REF = (fwd = 99.16322764791194, eu = 8.578451710277706, barrier = 7.561814171846508, asian_float = 4.86386704425551)

@testset "Normal Inverse Gaussian" begin
    theta1, k1 = 0.01, 0.03
    model = NormalInverseGaussianProcess(sigma, theta1, k1, Underlying(S0, d))

    @testset "standard" begin
        check_prices(model, rfCurve, mc, standard_payoffs(), NIG_REF; toll = 0.8)
    end

    @testset "antithetic" begin
        check_prices(model, rfCurve, mc_anti, standard_payoffs(), NIG_REF; toll = 0.8)
    end

    @testset "parameter validation" begin
        @test_throws AssertionError NormalInverseGaussianProcess(-sigma, theta1, k1, Underlying(S0, d))
        @test_throws AssertionError NormalInverseGaussianProcess(sigma, theta1, -k1, Underlying(S0, d))
        @test_throws AssertionError NormalInverseGaussianProcess(sigma, 10000.0, k1, Underlying(S0, d))
    end
end
