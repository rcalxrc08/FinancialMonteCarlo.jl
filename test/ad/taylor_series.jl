include("../testutils.jl")
using TaylorSeries

@testset "TaylorSeries" begin
    model = BlackScholesProcess(taylor_expand(identity, sigma, order = 3), Underlying(S0, d))
    check_prices(model, rfCurve, mc, standard_payoffs(), BS_CALL_REF; toll = 0.8)
    check_prices(model, rfCurve, mc_anti, standard_payoffs(), BS_CALL_REF; toll = 0.6)
end
