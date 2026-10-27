include("../testutils.jl")

# With (almost) all the weight on one component, both mixtures collapse to Black-Scholes.
@testset "LogNormal mixtures" begin
    eta = [0.2, 0.2]
    lam = Float64[0.999999999]
    und = Underlying(S0, d)

    @testset "$name" for (name, model, build) in [
        ("LogNormalMixture", LogNormalMixture(eta, lam, und), (η, λ) -> LogNormalMixture(η, λ, und)),
        ("ShiftedLogNormalMixture", ShiftedLogNormalMixture(eta, lam, 0.0, und), (η, λ) -> ShiftedLogNormalMixture(η, λ, 0.0, und)),
    ]
        @testset "parameters" begin
            params = FinancialMonteCarlo.get_parameters(model)
            @test FinancialMonteCarlo.get_parameters(FinancialMonteCarlo.set_parameters!(model, params)) == params
            @test occursin(name, capture_stdout(() -> display(model)))
        end

        @testset "prices" begin
            check_prices(model, rfCurve, mc, standard_payoffs(), BS_CALL_REF; toll = 0.8)
        end

        @testset "parameter validation" begin
            @test_throws AssertionError build(-0.2 * ones(3), 0.1 * ones(2))
            @test_throws AssertionError build(0.2 * ones(3), -0.1 * ones(2))
            @test_throws AssertionError build(0.2 * ones(3), ones(2))
            @test_throws AssertionError build(0.2 * ones(3), 0.2 * ones(3))
        end
    end
end
