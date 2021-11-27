using Test, FinancialMonteCarlo, DualNumbers;
@show "Black Scholes Model"
S0 = 100.0;
K = 100.0;
T = 1.0;
d = 0.01;
D = 90.0;

Nsim = 10000;
Nstep = 30;
sigma = 0.2;
mc = MonteCarloConfiguration(Nsim, Nstep);
dr__ = 1e-7;
mc1 = MonteCarloConfiguration(Nsim, Nstep, FinancialMonteCarlo.AntitheticMC());
toll = 0.8
der_toll = 1e-2;

rfCurve = ZeroRate(dual(0.02, 1.0));
rfCurve1 = ZeroRate(0.02);

FwdData = Forward(T)
EUData = EuropeanOption(T, K)
AMData = AmericanOption(T, K)
BarrierData = BarrierOptionDownOut(T, K, D)
AsianFloatingStrikeData = AsianFloatingStrikeOption(T)
AsianFixedStrikeData = AsianFixedStrikeOption(T, K)
Model = BlackScholesProcess(sigma, Underlying(S0, d));

@show EuPrice = pricer(Model, rfCurve, mc, EUData);
@show EuPrice_der = rho(Model, rfCurve1, mc, EUData, dr__);

@show AmPrice = pricer(Model, rfCurve, mc, AMData);
@show AmPrice_der = rho(Model, rfCurve1, mc, AMData, dr__);

@test abs(EuPrice - 8.43005524824866) < toll
@test abs(EuPrice.epsilon - EuPrice_der) < der_toll

@test abs(AmPrice - 8.450489415187354) < toll
@test abs(AmPrice.epsilon - AmPrice_der) < der_toll

@show EuPrice = pricer(Model, rfCurve, mc1, EUData);
@show EuPrice_der = rho(Model, rfCurve1, mc1, EUData, dr__);
@show EuPrice_der_vec = rho(Model, rfCurve1, mc1, [EUData], dr__)[1];
tollanti = 0.6;

@test abs(EuPrice - 8.43005524824866) < tollanti
@test abs(EuPrice.epsilon - EuPrice_der) < der_toll
@test abs(EuPrice.epsilon - EuPrice_der_vec) < der_toll
