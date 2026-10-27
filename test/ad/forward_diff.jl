include("../testutils.jl")
using ForwardDiff

@testset "ForwardDiff" begin
    @testset "Dual-valued volatility" begin
        model = BlackScholesProcess(ForwardDiff.Dual{Float64}(sigma, 1.0), Underlying(S0, d))
        # American exercise on dual paths is not checked against the reference with standard MC
        check_prices(model, rfCurve, mc, standard_payoffs(), Base.structdiff(BS_CALL_REF, NamedTuple{(:am,)}); toll = 0.8)
        check_prices(model, rfCurve, mc_anti, standard_payoffs(), BS_CALL_REF; toll = 0.6)
    end

    @testset "gradient of a European call" begin
        f(x) = pricer(BlackScholesProcess(x[1], Underlying(x[2], d)), rfCurve, mc, EuropeanOption(T, K))
        vega, delta_ = ForwardDiff.gradient(f, [sigma, S0])
        @test vega > 0
        @test 0 < delta_ < 1
    end
end
