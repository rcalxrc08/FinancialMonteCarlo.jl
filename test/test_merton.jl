using Test, FinancialMonteCarlo;
@show "MertonProcess"
S0 = 100.0;
K = 100.0;
r = 0.02;
T = 1.0;
d = 0.01;
D = 90.0;

Nsim = 10000;
Nstep = 30;
sigma = 0.2;
lam = 5.0;
mu1 = 0.03;
sigma1 = 0.02;
mc = MonteCarloConfiguration(Nsim, Nstep);
mc1 = MonteCarloConfiguration(Nsim, Nstep, FinancialMonteCarlo.AntitheticMC());
toll = 0.8;

rfCurve = ZeroRate(r);

FwdData = Forward(T)
EUData = EuropeanOption(T, K)
AMData = AmericanOption(T, K)
BarrierData = BarrierOptionDownOut(T, K, D)
AsianData = AsianFloatingStrikeOption(T)
Model = MertonProcess(sigma, lam, mu1, sigma1, Underlying(S0, d));

@show FwdPrice = pricer(Model, rfCurve, mc, FwdData);
@show EuPrice = pricer(Model, rfCurve, mc, EUData);
@show AmPrice = pricer(Model, rfCurve, mc, AMData);
@show BarrierPrice = pricer(Model, rfCurve, mc, BarrierData);
@show AsianPrice = pricer(Model, rfCurve, mc, AsianData);

@test abs(FwdPrice - 99.1188767166039) < toll
@test abs(EuPrice - 9.084327245917533) < toll
@test abs(BarrierPrice - 7.880881290426765) < toll
@test abs(AsianPrice - 5.129020349580892) < toll

@show FwdPrice = pricer(Model, rfCurve, mc1, FwdData);
@show EuPrice = pricer(Model, rfCurve, mc1, EUData);
@show AmPrice = pricer(Model, rfCurve, mc1, AMData);
@show BarrierPrice = pricer(Model, rfCurve, mc1, BarrierData);
@show AsianPrice = pricer(Model, rfCurve, mc1, AsianData);
tollanti = 0.6;
@test abs(FwdPrice - 99.1188767166039) < tollanti
@test abs(EuPrice - 9.084327245917533) < tollanti
@test abs(BarrierPrice - 7.880881290426765) < tollanti
@test abs(AsianPrice - 5.129020349580892) < tollanti

@show "Test Merton Parameters"

@test_throws(ErrorException, MertonProcess(-sigma, lam, mu1, sigma1, Underlying(S0, d)))
@test_throws(ErrorException, MertonProcess(sigma, lam, mu1, -sigma1, Underlying(S0, d)))
@test_throws(ErrorException, MertonProcess(sigma, -lam, mu1, sigma1, Underlying(S0, d)))
