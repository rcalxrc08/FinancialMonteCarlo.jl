include("../testutils.jl")
using Distributed, VectorizedRNG

@testset "Parallel modes and RNGs" begin
    model = BlackScholesProcess(sigma, Underlying(S0, d))
    eu = EuropeanOption(T, K)

    @testset "$name" for (name, mode) in [
        ("MultiThreading", FinancialMonteCarlo.MultiThreading()),
        ("MultiProcess", FinancialMonteCarlo.MultiProcess()),
        ("VectorizedRNG", FinancialMonteCarlo.SerialMode(10, local_rng())),
    ]
        mc_par = MonteCarloConfiguration(Nsim, Nstep, mode)
        @test abs(pricer(model, rfCurve, mc_par, eu) - BS_CALL_REF.eu) < 1.0
        prices = pricer(model, rfCurve, mc_par, [eu, eu])
        @test prices[1] ≈ prices[2]
    end
end
