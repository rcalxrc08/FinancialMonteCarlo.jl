using BenchmarkTools
using FinancialMonteCarlo, CUDA
CUDA.allowscalar(false)
S0 = 100.0;
K = 100.0;
r = 0.02;
T = 1.0;
d = 0.01;
D = 90.0;

Nsim = 10000;
Nstep = 30;
sigma = 0.2;
p = 0.3;
lam = 5.0;
lamp = 30.0;
lamm = 20.0;
mc = MonteCarloConfiguration(Nsim, Nstep);

mc_2 = MonteCarloConfiguration(Nsim, Nstep, FinancialMonteCarlo.CudaMode());
toll = 0.8;

rfCurve = ZeroRate(r);

FwdData = Forward(T)
EUData = EuropeanOption(T, K)
AMData = AmericanOption(T, K)
BarrierData = BarrierOptionDownOut(T, K, D)
AsianFloatingStrikeData = AsianFloatingStrikeOption(T)
AsianFixedStrikeData = AsianFixedStrikeOption(T, K)
Model = KouProcess(sigma, lam, p, lamp, lamm, Underlying(S0, d));

@show "CUDA_1 fwd"

@show "STD fwd"
@btime FwdPrice = pricer(Model, rfCurve, mc, FwdData);
@show "CUDA_1 fwd"

@show "CUDA_2 fwd"
@btime FwdPrice = pricer(Model, rfCurve, mc_2, FwdData);
@show "std eu"
@btime EuPrice = pricer(Model, rfCurve, mc, EUData);
@show "CUDA_1 eu"

@show "CUDA_2 eu"
@btime EuPrice = pricer(Model, rfCurve, mc_2, EUData);
@show "std am"
@btime AmPrice = pricer(Model, rfCurve, mc, AMData);
@show "CUDA_1 am"

@show "CUDA_2 am"
@btime AmPrice = pricer(Model, rfCurve, mc_2, AMData);
