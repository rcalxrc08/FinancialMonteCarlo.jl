using Test, FinancialMonteCarlo, DualNumbers;
@show "Black Scholes Model"
S0 = 100.0;
K = 100.0;
r = 0.02;
T = 1.0;
d = 0.01;
D = 90.0;

underlying_name = "ENJ"
underlying_name_2 = "ABPL"

underlying_ = underlying_name * "_" * underlying_name_2 * "_TESL"

Nsim = 10000;
Nstep = 30;
sigma = 0.2;
mc = MonteCarloConfiguration(Nsim, Nstep);
mc1 = MonteCarloConfiguration(Nsim, Nstep, FinancialMonteCarlo.AntitheticMC());
toll = 0.8

rfCurve = ZeroRate(r);

FwdData = Forward(T)
EUData = EuropeanOptionND(T, K)
AMData = AmericanOption(T, K)
BarrierData = BarrierOptionDownOut(T, K, D)
AsianFloatingStrikeData = AsianFloatingStrikeOption(T)
AsianFixedStrikeData = AsianFixedStrikeOption(T, K)
Model_enj = BlackScholesProcess(sigma, Underlying(S0, d));
Model_abpl = BlackScholesProcess(sigma, Underlying(S0, d));
Model_tesl = BlackScholesProcess(sigma, Underlying(S0, d));
rho_1 = [1.0 0.0 0.1; 0.0 1.0 0.0; 0.1 0.0 1.0];
Model = GaussianCopulaNVariateProcess(rho_1, Model_enj, Model_abpl, Model_tesl)

display(Model)

mktdataset = underlying_ → Model
mktdataset_1 = underlying_name → Model_enj

portfolio_ = [EUData];
portfolio = underlying_ → EUData
portfolio_1 = underlying_name → AMData
price_mkt = FinancialMonteCarlo.delta(mktdataset, rfCurve, mc, portfolio, underlying_name, 1e-7)
price_mkt_1 = FinancialMonteCarlo.delta(mktdataset_1, rfCurve, mc, portfolio_1, underlying_name, 1e-7)

Model_enj_dual = BlackScholesProcess(sigma, Underlying(dual(S0, 1.0), d));
Model_dual = GaussianCopulaNVariateProcess(rho_1, Model_enj_dual, Model_abpl, Model_tesl)

mktdataset_dual = underlying_ → Model_dual
mktdataset_1_dual = underlying_name → Model_enj_dual
price_mkt_dual = pricer(mktdataset_dual, rfCurve, mc, portfolio)
price_mkt_1_dual = pricer(mktdataset_1_dual, rfCurve, mc, portfolio_1)

@test(abs(price_mkt_dual.epsilon - price_mkt) < 1e-4)
@test(abs(price_mkt_1_dual.epsilon - price_mkt_1) < 1e-4)
@test(abs(price_mkt_1_dual.epsilon - FinancialMonteCarlo.delta(Model_enj, rfCurve, mc, AMData, 1e-7)) < 1e-4)
@test(abs(price_mkt_1_dual.epsilon - FinancialMonteCarlo.delta(Model_enj, rfCurve, mc, AMData * 1.0, 1e-7)) < 1e-4)
@test(abs(price_mkt_1_dual.epsilon - FinancialMonteCarlo.delta(Model_enj, rfCurve, mc, [AMData], 1e-7)[1]) < 1e-4)

@test_throws(AssertionError, FinancialMonteCarlo.delta(mktdataset_1, rfCurve, mc, portfolio_1, underlying_, 1e-7))
@test(FinancialMonteCarlo.delta(mktdataset_1, rfCurve, mc, portfolio_1, "ciao", 1e-7) == 0.0)
