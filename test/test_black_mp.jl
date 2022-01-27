using FinancialMonteCarlo, Distributed
S0 = 100.0;
K = 100.0;
r = 0.02;
T = 1.0;
d = 0.01;
D = 90.0;
Nsim = 10000;
Nstep = 30;
sigma = 0.2;
mc = MonteCarloConfiguration(Nsim, Nstep, FinancialMonteCarlo.MultiProcess());
toll = 1.0;
rfCurve = ZeroRate(r);
FwdData = Forward(T)
EUData = EuropeanOption(T, K)
AMData = AmericanOption(T, K)
BarrierData = BarrierOptionDownOut(T, K, D)
AsianFloatingStrikeData = AsianFloatingStrikeOption(T)
AsianFixedStrikeData = AsianFixedStrikeOption(T, K)
Model = BlackScholesProcess(sigma, Underlying(S0, d));
EuPrice = pricer(Model, rfCurve, mc, EUData);

@test abs(EuPrice - 8.43005524824866) < toll