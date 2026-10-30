include("../testutils.jl")
using DifferentialEquations, Random

# Pricing straight from an SDE definition through the DifferentialEquations.jl extension.
@testset "DifferentialEquations models" begin
    @testset "Black-Scholes SDE" begin
        drift(u, p, t) = (r - d) * u
        diffusion(u, p, t) = sigma * u
        prob = EnsembleProblem(SDEProblem(drift, diffusion, S0, (0.0, T)))
        model = FinancialMonteCarlo.MonteCarloDiffEqModel(prob, Underlying(S0, d))
        check_prices(model, rfCurve, mc_small, standard_payoffs())
    end

    @testset "Kou jump diffusion (log price)" begin
        p1, lam, lamp, lamm = 0.3, 5.0, 30.0, 20.0
        rng = mc.parallelMode.rng
        drift(u, p, t) = r - d - sigma^2 / 2.0 - lam * (p1 / (lamp - 1) - (1 - p1) / (lamm + 1))
        diffusion(u, p, t) = sigma
        rate(u, p, t) = lam * T
        affect!(integrator) = (integrator.u += rand(rng) < p1 ? randexp(rng) / lamp : -randexp(rng) / lamm)
        jump_prob = JumpProblem(SDEProblem(drift, diffusion, 0.0, (0.0, T)), Direct(), ConstantRateJump(rate, affect!), rng = rng)
        model = FinancialMonteCarlo.MonteCarloDiffEqModel(EnsembleProblem(jump_prob), x -> S0 * exp(x), Underlying(S0, d))
        payoffs = standard_payoffs()
        check_prices(model, rfCurve, mc, (fwd = payoffs.fwd,), (fwd = 98.8436678850961,); toll = 0.8)
        # the jumps draw from rng, which set_seed! must reset here too
        mc_small_jumps = MonteCarloConfiguration(Nsim_small, Nstep_small, FinancialMonteCarlo.StandardMC(), FinancialMonteCarlo.SerialMode(0, rng))
        check_prices(model, rfCurve, mc_small_jumps, Base.structdiff(payoffs, NamedTuple{(:fwd,)}))
    end
end
