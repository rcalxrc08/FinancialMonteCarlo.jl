include("../testutils.jl")
using TaylorSeries

@testset "TaylorSeries" begin
    model = BlackScholesProcess(taylor_expand(identity, sigma, order = 3))
    check_prices(model, und, rfCurve, mc, standard_payoffs(), BS_CALL_REF; toll = 0.8)
    check_prices(model, und, rfCurve, mc_small_anti, standard_payoffs())
    # max(::Taylor1, ::Float64) is the fallback of Base, which returns either argument unconverted: the payoffs compare
    # with a zero of the type of the paths, so they stay inferred (#53)
    x = taylor_expand(identity, 1.0, order = 2)
    paths = [fill(x * S0, Nsim_small, Nstep_small + 1) for _ = 1:2]
    @test eltype(@inferred payoff(paths, EuropeanOptionND(T, 3 * K), rfCurve, mc_small)) == Taylor1{Float64}
end
