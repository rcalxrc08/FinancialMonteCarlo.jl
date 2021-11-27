using BenchmarkTools
using FinancialMonteCarlo, CUDA
CUDA.allowscalar(false)
S0 = 100.0f0;
K = 100.0f0;
r = 0.02f0;
T = 1.0f0;
d = 0.01f0;
D = 90.0f0;

Nsim = 10000;
Nstep = 30;
sigma = 0.2f0;
mc = MonteCarloConfiguration(Nsim, Nstep);
mc_2 = MonteCarloConfiguration(Nsim, Nstep, FinancialMonteCarlo.CudaMode());
toll = 1e-3;

rfCurve = ZeroRate(r);

FwdData = Forward(T)
EUData = EuropeanOption(T, K)
AMData = AmericanOption(T, K)
BarrierData = BarrierOptionDownOut(T, K, D)
AsianFloatingStrikeData = AsianFloatingStrikeOption(T)
AsianFixedStrikeData = AsianFixedStrikeOption(T, K)
Model = BlackScholesProcess(sigma, Underlying(S0, d));

@show "CUDA_2 fwd"
@btime FwdPrice = pricer(Model, rfCurve, mc_2, FwdData);
@show "CUDA_2 eu"
@btime EuPrice = pricer(Model, rfCurve, mc_2, EUData);
