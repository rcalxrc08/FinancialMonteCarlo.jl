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
const mc = MonteCarloConfiguration(Nsim, Nstep);
toll = 1e-3;

const rfCurve = ZeroRate(r);

const FwdData = Forward(T)
const EUData = EuropeanOption(T, K)
const AMData = AmericanOption(T, K)
const BarrierData = BarrierOptionDownOut(T, K, D)
const AsianFloatingStrikeData = AsianFloatingStrikeOption(T)
const AsianFixedStrikeData = AsianFixedStrikeOption(T, K)
const Model = BlackScholesProcess(sigma, Underlying(S0, d));

@btime FwdPrice = pricer(Model, rfCurve, mc, FwdData);
@btime EuPrice = pricer(Model, rfCurve, mc, EUData);
@btime AmPrice = pricer(Model, rfCurve, mc, AMData);
@btime BarrierPrice = pricer(Model, rfCurve, mc, BarrierData);
@btime AsianPrice1 = pricer(Model, rfCurve, mc, AsianFloatingStrikeData);
@btime AsianPrice1 = pricer(Model, rfCurve, mc, AsianFixedStrikeData);

const optionDatas = [FwdData, EUData, AMData, BarrierData, AsianFloatingStrikeData, AsianFixedStrikeData]

@btime (FwdPrice, EuPrice, AMPrice, BarrierPrice, AsianPrice1, AsianPrice2) = pricer(Model, rfCurve, mc, optionDatas)
