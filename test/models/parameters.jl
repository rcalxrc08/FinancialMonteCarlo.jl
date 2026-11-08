include("../testutils.jl")
using ForwardDiff

# get_parameters / set_parameters! round trip for every model; both must be inferred (callers such as
# calibration loops rely on it), with Float64 and with dual parameters.
@testset "get_parameters and set_parameters!" begin
    und = Underlying(S0, d)
    models = [BlackScholesProcess(sigma, und), HestonProcess(0.2, 0.2, 0.1, 1.5, -0.5, 0.04, und), KouProcess(0.2, 3.0, 0.3, 9.0, 5.0, und), MertonProcess(0.2, 3.0, -0.05, 0.1, und), VarianceGammaProcess(0.2, -0.03, 0.16, und), NormalInverseGaussianProcess(0.2, -0.03, 0.16, und), LogNormalMixture([0.2, 0.3, 0.25], [0.3, 0.4], und), ShiftedLogNormalMixture([0.2, 0.3, 0.25], [0.3, 0.4], 0.1, und), BrownianMotion(0.2, 0.01)]
    @testset "$(nameof(typeof(model)))" for model in models
        params = @inferred FinancialMonteCarlo.get_parameters(model)
        new_model = @inferred FinancialMonteCarlo.set_parameters!(model, params)
        @test FinancialMonteCarlo.get_parameters(new_model) == params
        model isa FinancialMonteCarlo.AbstractMonteCarloProcess && @test new_model.underlying === model.underlying
        dual_params = ForwardDiff.Dual{Nothing}.(params, 1.0)
        dual_model = @inferred FinancialMonteCarlo.set_parameters!(model, dual_params)
        @test (@inferred FinancialMonteCarlo.get_parameters(dual_model)) == dual_params
        @test_throws AssertionError FinancialMonteCarlo.set_parameters!(model, [params; params])
    end
end

@testset "get_parameters of payoffs and copulas" begin
    @test (@inferred FinancialMonteCarlo.get_parameters(EuropeanOption(T, K))) == [T, K, 1.0]
    @test (@inferred FinancialMonteCarlo.get_parameters(EuropeanOption(Float32(T), Float32(K)))) isa Vector{Float32}
    @test (@inferred FinancialMonteCarlo.get_parameters(BermudanOption([T / 2, T], K))) == [T, T / 2, T, K, 1.0]
    @test isempty(@inferred FinancialMonteCarlo.get_parameters(Spot()))
    bs = BlackScholesProcess(sigma, Underlying(S0, d))
    @test (@inferred FinancialMonteCarlo.get_parameters(GaussianCopulaNVariateProcess(bs, bs, 0.4))) == [1.0, 0.4, 0.4, 1.0]
end
