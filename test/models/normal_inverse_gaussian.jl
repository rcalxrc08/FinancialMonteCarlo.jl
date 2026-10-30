include("../testutils.jl")
using ForwardDiff

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

    @testset "drift and subordinator law with explicit parameters" begin
        dt = T / Nstep
        @test FinancialMonteCarlo.compute_drift(model) == FinancialMonteCarlo.compute_drift(model, sigma, theta1, k1)
        @test FinancialMonteCarlo.compute_subordinator_law(model, dt) == FinancialMonteCarlo.compute_subordinator_law(model, k1, dt)
        @test FinancialMonteCarlo.compute_drift(model, Float32(sigma), Float32(theta1), Float32(k1)) isa Float32
        @test FinancialMonteCarlo.compute_subordinator_law(model, Float32(k1), Float32(dt)) isa FinancialMonteCarlo.InverseGaussian{Float32}
        model_dual = NormalInverseGaussianProcess(sigma, theta1, ForwardDiff.Dual(k1, 1.0), Underlying(S0, d))
        @test FinancialMonteCarlo.compute_subordinator_law(model_dual, dt) isa FinancialMonteCarlo.InverseGaussian{<:ForwardDiff.Dual}
        @test FinancialMonteCarlo.compute_subordinator_law(model_dual, ForwardDiff.value(model_dual.κ), dt) isa FinancialMonteCarlo.InverseGaussian{Float64}
    end

    @testset "parameter validation" begin
        @test_throws AssertionError NormalInverseGaussianProcess(-sigma, theta1, k1, Underlying(S0, d))
        @test_throws AssertionError NormalInverseGaussianProcess(sigma, theta1, -k1, Underlying(S0, d))
        @test_throws AssertionError NormalInverseGaussianProcess(sigma, 10000.0, k1, Underlying(S0, d))
    end
end
