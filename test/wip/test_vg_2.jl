using Test, FinancialMonteCarlo, RandomNumbers, DualNumbers;
@show "VarianceGammaProcess"
S0 = dual(100.0, 1.0);
K = 105.0;
r = 0.02;
T = 1.0;
d = 0.01;
D = 90.0;

Nsim = 100000;
Nstep = 300;
sigma = 0.2;
theta1 = 0.01;
k1 = 1.03;
#mc = MonteCarloConfiguration(Nsim, Nstep, FinancialMonteCarlo.SerialMode(1, RandomNumbers.Xorshifts.Xoroshiro128Plus()));
mc = MonteCarloConfiguration(Nsim, Nstep)#, FinancialMonteCarlo.SerialMode(1, RandomNumbers.Xorshifts.Xoroshiro128Plus()));
mc1 = MonteCarloConfiguration(Nsim, Nstep, FinancialMonteCarlo.AntitheticMC());
toll = 0.8;

rfCurve = ZeroRate(r);

FwdData = Forward(T)
EUData = EuropeanOption(T, K)
AMData = AmericanOption(T, K)
BarrierData = BarrierOptionDownOut(T, K, D)
AsianData = AsianFloatingStrikeOption(T)
Model = NormalInverseGaussianProcess(sigma, theta1, k1, Underlying(S0, d));
@show EuPrice = pricer(Model, rfCurve, mc, EUData);
@show EuPrice = pricer(Model, rfCurve, mc1, EUData);
