using BenchmarkTools
using FinancialMonteCarlo

S0 = 100.0;
K = 100.0;
r = 0.02;
T = 1.0;
d = 0.01;
D = 90.0;

Nsim = 10000;
Nstep = 30;
sigma = 0.2;
mc = MonteCarloConfiguration(Nsim, Nstep, FinancialMonteCarlo.AntitheticMC());
toll = 1e-3;

rfCurve = ZeroRate(r);

FwdData = Forward(T)
EUData = EuropeanOption(T, K)
AMData = AmericanOption(T, K)
BarrierData = BarrierOptionDownOut(T, K, D)
AsianFloatingStrikeData = AsianFloatingStrikeOption(T)
AsianFixedStrikeData = AsianFixedStrikeOption(T, K)
Model = BlackScholesProcess(sigma, Underlying(S0, d));

@btime FwdPrice = pricer(Model, rfCurve, mc, FwdData);
@btime EuPrice = pricer(Model, rfCurve, mc, EUData);
@btime AmPrice = pricer(Model, rfCurve, mc, AMData);
@btime BarrierPrice = pricer(Model, rfCurve, mc, BarrierData);
@btime AsianPrice1 = pricer(Model, rfCurve, mc, AsianFloatingStrikeData);
@btime AsianPrice1 = pricer(Model, rfCurve, mc, AsianFixedStrikeData);

optionDatas = [FwdData, EUData, AMData, BarrierData, AsianFloatingStrikeData, AsianFixedStrikeData]

@btime (FwdPrice, EuPrice, AMPrice, BarrierPrice, AsianPrice1, AsianPrice2) = pricer(Model, rfCurve, mc, optionDatas)
