using Test, FinancialMonteCarlo, ForwardDiff;
@show "Black Scholes Model"
S0 = 100.0;
K = 100.0;
r = 0.02;
T = 1.0;
d = 0.01;
D = 90.0;

Nsim = 10000;
Nstep = 30;
sigma = ForwardDiff.Dual{Float64}(0.2, 1.0)
mc = MonteCarloConfiguration(Nsim, Nstep);
mc1 = MonteCarloConfiguration(Nsim, Nstep, FinancialMonteCarlo.AntitheticMC());
toll = 0.8

rfCurve = ZeroRate(r);

FwdData = Forward(T)
EUData = EuropeanOption(T, K)
AMData = AmericanOption(T, K)
BarrierData = BarrierOptionDownOut(T, K, D)
AsianFloatingStrikeData = AsianFloatingStrikeOption(T)
AsianFixedStrikeData = AsianFixedStrikeOption(T, K)
Model = BlackScholesProcess(sigma, Underlying(S0, d));

@show FwdPrice = pricer(Model, rfCurve, mc, FwdData);
@show EuPrice = pricer(Model, rfCurve, mc, EUData);
@show AmPrice = pricer(Model, rfCurve, mc, AMData);
@show BarrierPrice = pricer(Model, rfCurve, mc, BarrierData);
@show AsianPrice1 = pricer(Model, rfCurve, mc, AsianFloatingStrikeData);
@show AsianPrice2 = pricer(Model, rfCurve, mc, AsianFixedStrikeData);

@test abs(FwdPrice - 99.1078451563562) < toll
@test abs(EuPrice - 8.43005524824866) < toll
@test abs(BarrierPrice - 7.5008664470880735) < toll
@test abs(AsianPrice1 - 4.774451704549382) < toll

@show FwdPrice = pricer(Model, rfCurve, mc1, FwdData);
@show EuPrice = pricer(Model, rfCurve, mc1, EUData);
@show AmPrice = pricer(Model, rfCurve, mc1, AMData);
@show BarrierPrice = pricer(Model, rfCurve, mc1, BarrierData);
@show AsianPrice1 = pricer(Model, rfCurve, mc1, AsianFloatingStrikeData);
@show AsianPrice2 = pricer(Model, rfCurve, mc1, AsianFixedStrikeData);
tollanti = 0.6;
@test abs(FwdPrice - 99.1078451563562) < tollanti
@test abs(EuPrice - 8.43005524824866) < tollanti
@test abs(AmPrice - 8.450489415187354) < tollanti
@test abs(BarrierPrice - 7.5008664470880735) < tollanti
@test abs(AsianPrice1 - 4.774451704549382) < tollanti
