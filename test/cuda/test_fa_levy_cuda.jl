using Test, FinancialMonteCarlo, CUDA;
const FMC = FinancialMonteCarlo

S0 = 100.0;
K = 100.0;
r = 0.02;
T = 1.0;
d = 0.01;

Nsim = 100000;
Nstep = 30;
toll = 0.35;
rfCurve = ZeroRate(r);
EUData = EuropeanOption(T, K);
und = Underlying(S0, d);
models = [KouProcess(0.2, 3.0, 0.5, 25.0, 30.0, und), MertonProcess(0.2, 3.0, 0.03, 0.02, und)];

@testset "$(nameof(typeof(model)))" for model in models
    @testset "$(nameof(typeof(method)))" for method in (FMC.StandardMC(), FMC.AntitheticMC())
        mc = MonteCarloConfiguration(Nsim, Nstep, method, FMC.CudaMode(1, 2))
        mc_cpu = MonteCarloConfiguration(Nsim, Nstep, method, FMC.SerialMode(3))
        # The generic method also works on a CuMatrix, launching a kernel for each jump.
        sim_method = which(FMC.simulate!, Tuple{CuMatrix{Float64}, typeof(model), typeof(rfCurve), typeof(mc), Float64})
        @test sim_method.module === Base.get_extension(FMC, :FinancialMonteCarloCUDAExt)
        S = FMC.simulate(model, rfCurve, mc, T)
        @test S isa CuMatrix{Float64}
        @test abs(pricer(model, rfCurve, mc, EUData) - pricer(model, rfCurve, mc_cpu, EUData)) < toll
    end
end
