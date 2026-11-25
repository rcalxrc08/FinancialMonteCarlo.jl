include("../testutils.jl")
using FinancialToolbox: blsdownout, blsdownin, blsupout, blsupin
using ForwardDiff

# The single-barrier options under Black-Scholes against the closed forms of FinancialToolbox (#63). Those are
# continuously monitored, the Monte Carlo price only on the Nstep + 1 dates of the grid, which biases it towards the
# knock-out side by O(√Δt). Broadie, Glasserman and Kou: shifting the barrier of the closed form away from the spot by
# exp(β σ √Δt), β = -ζ(1/2)/√(2π) ≈ 0.5826, gives the price of the discrete one to o(√Δt).
@testset "Barrier options: closed forms" begin
    Nsim_bar, Nstep_bar = 100_000, 50
    mc_bar = MonteCarloConfiguration(Nsim_bar, Nstep_bar; sampling = FinancialMonteCarlo.AntitheticMC(), seed = 1)
    model = BlackScholesProcess(sigma)
    shift = exp(0.5826 * sigma * sqrt(T / Nstep_bar))
    H_down, H_up = D, 1.15 * S0
    # option, closed form, barrier, factor of the discrete barrier, knock-out
    cases = [(BarrierOptionDownOut, blsdownout, H_down, 1 / shift, true), (BarrierOptionDownIn, blsdownin, H_down, 1 / shift, false), (BarrierOptionUpOut, blsupout, H_up, shift, true), (BarrierOptionUpIn, blsupin, H_up, shift, false)]
    # strikes below the down barrier, between the barriers and above the up barrier
    @testset "$(nameof(Option)), K = $strike, $(isCall ? "call" : "put")" for (Option, closed_form, H, factor, knock_out) in cases, strike in (0.85 * S0, S0, 1.2 * S0), isCall in (true, false)
        payoff_ = Option(T, strike, H, isCall)
        price = pricer(model, und, rfCurve, mc_bar, payoff_)
        std_error = sqrt(variance(model, und, rfCurve, mc_bar, payoff_) / Nsim_bar)
        discrete = closed_form(S0, strike, r, T, sigma, d, H * factor, isCall)
        continuous = closed_form(S0, strike, r, T, sigma, d, H, isCall)
        # the correction is o(√Δt), not exact: 2% of the price on top of the statistical error
        @test abs(price - discrete) <= 4 * std_error + discrete / 50
        # monitoring on a grid misses crossings: knock-out options are worth more, knock-in options less
        @test knock_out ? price >= continuous - 4 * std_error : price <= continuous + 4 * std_error
    end

    # Pathwise derivatives (AD) differentiate the payoff path by path, so they miss the paths that cross the barrier
    # when the spot or the volatility moves: the delta of a down-and-out call at the money is two thirds of the closed form (#74)
    @testset "pathwise delta" begin
        price_of_spot(S) = pricer(BlackScholesProcess(sigma), Underlying(S, d), rfCurve, mc_bar, BarrierOptionDownOut(T, K, H_down))
        delta_closed_form = ForwardDiff.derivative(S -> blsdownout(S, K, r, T, sigma, d, H_down / shift), S0)
        @test_broken abs(ForwardDiff.derivative(price_of_spot, S0) - delta_closed_form) < 0.05
    end
end
