using Test, FinancialMonteCarlo;
@show "Black Scholes Model"
S0 = 100.0;
K = 100.0;
r = 0.02;
T = 1.0;
d = 0.01;
D = 90.0;

Nsim = 10000;
Nstep = 30;
sigma = 0.2;
mc = MonteCarloConfiguration(Nsim, Nstep);
mc1 = MonteCarloConfiguration(Nsim, Nstep, FinancialMonteCarlo.AntitheticMC());
toll = 0.8

rfCurve = ZeroRate(r);

FwdData = Forward(T)
EUData = EuropeanOption(T, K)
EUBin = BinaryEuropeanOption(T, K)
AMData = AmericanOption(T, K)
AmBin = BinaryEuropeanOption(T, K)
BarrierData = BarrierOptionDownOut(T, K, D)
AsianFloatingStrikeData = AsianFloatingStrikeOption(T)
AsianFixedStrikeData = AsianFixedStrikeOption(T, K)
Model = BlackScholesProcess(sigma, Underlying(S0, d));

display(Model)
ST = FinancialMonteCarlo.distribution(Model, rfCurve, mc, T);
@show FwdPrice = exp(-r * T) * sum(ST) / length(ST);

@test abs(FwdPrice - 99.1078451563562) < toll
