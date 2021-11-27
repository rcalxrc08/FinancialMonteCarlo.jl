using BenchmarkTools
using FinancialMonteCarlo

S0 = 100.0;
K = 100.0;
r = [0.00, 0.02];
T = 1.0;
d = FinancialMonteCarlo.Curve([0.00, 0.02], T);
D = 90.0;

Nsim = 10000;
Nstep = 30;
sigma = 0.2;
p = 0.3;
lam = 5.0;
lamp = 30.0;
lamm = 20.0;
mc = MonteCarloConfiguration(Nsim, Nstep);
toll = 0.8;

rfCurve = FinancialMonteCarlo.ZeroRateCurve(r, T);

FwdData = Forward(T)
EUData = EuropeanOption(T, K)
AMData = AmericanOption(T, K)
BarrierData = BarrierOptionDownOut(T, K, D)
AsianFloatingStrikeData = AsianFloatingStrikeOption(T)
AsianFixedStrikeData = AsianFixedStrikeOption(T, K)
Model = KouProcess(sigma, lam, p, lamp, lamm, Underlying(S0, d));

@btime FwdPrice = pricer(Model, rfCurve, mc, FwdData);
@btime EuPrice = pricer(Model, rfCurve, mc, EUData);
@btime AmPrice = pricer(Model, rfCurve, mc, AMData);
@btime BarrierPrice = pricer(Model, rfCurve, mc, BarrierData);
@btime AsianPrice1 = pricer(Model, rfCurve, mc, AsianFloatingStrikeData);
@btime AsianPrice2 = pricer(Model, rfCurve, mc, AsianFixedStrikeData);

optionDatas = [FwdData, EUData, AMData, BarrierData, AsianFloatingStrikeData, AsianFixedStrikeData]

@btime (FwdPrice, EuPrice, AMPrice, BarrierPrice, AsianPrice1, AsianPrice2) = pricer(Model, rfCurve, mc, optionDatas)
