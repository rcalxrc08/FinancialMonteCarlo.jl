using Test, FinancialMonteCarlo, RandomNumbers;
@show "VarianceGammaProcess"
S0 = 100.0;
K = 100.0;
r = 0.02;
T = 1.0;
d = 0.01;
D = 90.0;

Nsim = 10000;
Nstep = 30;
sigma = 0.2;
theta1 = 0.01;
k1 = 0.03;
mc = MonteCarloConfiguration(Nsim, Nstep, FinancialMonteCarlo.SerialMode(1, RandomNumbers.Xorshifts.Xoroshiro128Plus()));
mc1 = MonteCarloConfiguration(Nsim, Nstep, FinancialMonteCarlo.AntitheticMC());
toll = 0.8;

rfCurve = ZeroRate(r);

FwdData = Forward(T)
EUData = EuropeanOption(T, K)
AMData = AmericanOption(T, K)
BarrierData = BarrierOptionDownOut(T, K, D)
AsianData = AsianFloatingStrikeOption(T)
Model = VarianceGammaProcess(sigma, theta1, k1, Underlying(S0, d));

@show FwdPrice = pricer(Model, rfCurve, mc, FwdData);
@show EuPrice = pricer(Model, rfCurve, mc, EUData);
@show AmPrice = pricer(Model, rfCurve, mc, AMData);
@show BarrierPrice = pricer(Model, rfCurve, mc, BarrierData);
@show AsianPrice = pricer(Model, rfCurve, mc, AsianData);

@test abs(FwdPrice - 99.09861035102614) < toll
@test abs(EuPrice - 8.368903690692187) < toll
@test abs(BarrierPrice - 7.469024475794258) < toll
@test abs(AsianPrice - 4.779663836736272) < toll

@show FwdPrice = pricer(Model, rfCurve, mc1, FwdData);
@show EuPrice = pricer(Model, rfCurve, mc1, EUData);
@show AmPrice = pricer(Model, rfCurve, mc1, AMData);
@show BarrierPrice = pricer(Model, rfCurve, mc1, BarrierData);
@show AsianPrice = pricer(Model, rfCurve, mc1, AsianData);

@test abs(FwdPrice - 99.09861035102614) < toll
@test abs(EuPrice - 8.368903690692187) < toll
@test abs(BarrierPrice - 7.469024475794258) < toll
@test abs(AsianPrice - 4.779663836736272) < toll
print(Model);
println("");
@show "Test Variance Gamma Parameters"

@test_throws(AssertionError, VarianceGammaProcess(-sigma, theta1, k1, Underlying(S0, d)))
@test_throws(AssertionError, VarianceGammaProcess(sigma, theta1, -k1, Underlying(S0, d)))
@test_throws(AssertionError, VarianceGammaProcess(sigma, 10000.0, k1, Underlying(S0, d)))
