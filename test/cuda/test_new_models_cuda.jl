using Test, FinancialMonteCarlo, CUDA;
const FMC = FinancialMonteCarlo

S0, K, r, T, d = 100.0, 100.0, 0.02, 1.0, 0.01
Nsim, Nstep, toll = 100_000, 30, 0.35
und = Underlying(S0, d)
payoffs = (EuropeanOption(T, K), BarrierOptionDownOut(T, K, 80.0), AsianFixedStrikeOption(T, K))
models = [BachelierProcess(20.0, und), LogNormalMixture((0.15, 0.3), (0.5,), und), ShiftedProcess(LogNormalMixture((0.15, 0.3), (0.5,), und), 0.2), ShiftedProcess(BlackScholesProcess(0.2, und), 0.3), ShiftedProcess(KouProcess(0.2, 3.0, 0.5, 25.0, 30.0, und), 0.2), ShiftedProcess(VarianceGammaProcess(0.2, -0.1, 0.2, und), 0.2)]

@testset "$(nameof(typeof(model)))" for model in models
    @testset "$(nameof(typeof(method))), $(typeof(rfCurve))" for method in (FMC.StandardMC(), FMC.AntitheticMC()), rfCurve in (ZeroRate(r), ZeroRate([0.01, 0.02], T))
        mc = MonteCarloConfiguration(Nsim, Nstep; sampling = method, mode = FMC.CudaMode(), seed = 1)
        mc_cpu = MonteCarloConfiguration(Nsim, Nstep; sampling = method, seed = 3)
        @test FMC.simulate(model, rfCurve, mc, T) isa CuMatrix{Float64}
        for payoff_ in payoffs
            @test abs(pricer(model, rfCurve, mc, payoff_) - pricer(model, rfCurve, mc_cpu, payoff_)) < toll
        end
    end
end
