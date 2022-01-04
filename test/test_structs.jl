using FinancialMonteCarlo, Test
@show "Test Structs Building"
S0 = 100.0;
K = 100.0;
Kneg = -100.0;
r = 0.02;
T = 1.0;
Tneg = -1.0;
d = 0.01;
D = 90.0;

Nsim = 10000;
Nstep = 30;
sigma = 0.2;
theta1 = 0.01;
k1 = 0.03;
sigma1 = 0.02;
underlying_name = "ENJ"
mc = MonteCarloConfiguration(Nsim, Nstep, 1);
mc2 = MonteCarloConfiguration(Nsim, Nstep, FinancialMonteCarlo.SerialMode(), 1);
mc1 = MonteCarloConfiguration(Nsim, Nstep, FinancialMonteCarlo.AntitheticMC());

toll = 0.8;
Model_enj = BlackScholesProcess(sigma, Underlying(S0, d));
Model_abpl = BlackScholesProcess(sigma, Underlying(S0, d));
FinancialMonteCarlo.Curve([0.0, 0.1], [0.0, 0.1]);
Model = GaussianCopulaNVariateProcess(Model_enj, Model_abpl, 0.0)
param_ = FinancialMonteCarlo.get_parameters(Model_abpl);

@test_throws(AssertionError, FinancialMonteCarlo.set_parameters!(Model_abpl, [1, 2, 3, 4]))

port_ = "eni_aap" → Model;
rfCurve = ZeroRate(r);
rfCurve2 = ZeroRate([0.00, 0.02], T);
@test_throws(AssertionError, port_ + port_)
@test_throws(AssertionError, underlying_name → Model)
@test_throws(AssertionError, Underlying(-S0, rfCurve2.r))
@test_throws(AssertionError, "c_i" → Forward(1.0))
@test_throws(AssertionError, Underlying(-S0))
@test_throws(AssertionError, EuropeanOptionND(-1.0, 1.0))
@test_throws(AssertionError, EuropeanOptionND(1.0, -1.0))
@test_throws(AssertionError, MonteCarloConfiguration(-Nsim, Nstep))
@test_throws(AssertionError, MonteCarloConfiguration(Nsim + 1, Nstep, FinancialMonteCarlo.AntitheticMC()));
@test_throws(AssertionError, MonteCarloConfiguration(Nsim, -Nstep))
#####################################################################
#Payoff Structs Test: Negative TTM
@test_throws(AssertionError, Forward(Tneg));
@test_throws(AssertionError, EuropeanOption(Tneg, K));
@test_throws(AssertionError, AmericanOption(Tneg, K));
@test_throws(AssertionError, BermudanOption([Tneg], K));
@test_throws(AssertionError, BinaryAmericanOption(Tneg, K));
@test_throws(AssertionError, BinaryEuropeanOption(Tneg, K));
@test_throws(AssertionError, BarrierOptionDownOut(Tneg, K, D));
@test_throws(AssertionError, BarrierOptionUpOut(Tneg, K, D));
@test_throws(AssertionError, BarrierOptionUpIn(Tneg, K, D));
@test_throws(AssertionError, BarrierOptionDownIn(Tneg, K, D));
@test_throws(AssertionError, AsianFloatingStrikeOption(Tneg));
@test_throws(AssertionError, AsianFixedStrikeOption(Tneg, K));
@test_throws(AssertionError, DoubleBarrierOption(Tneg, K, K * 10, D));

# Negative Strike and Barriers
@test_throws(AssertionError, EuropeanOption(T, Kneg));
@test_throws(AssertionError, AmericanOption(T, Kneg));
@test_throws(AssertionError, BermudanOption([T], Kneg));
@test_throws(AssertionError, BinaryEuropeanOption(T, Kneg));
@test_throws(AssertionError, BinaryAmericanOption(T, Kneg));
@test_throws(AssertionError, AsianFixedStrikeOption(T, Kneg));
@test_throws(AssertionError, BarrierOptionDownOut(T, Kneg, D));
@test_throws(AssertionError, BarrierOptionDownOut(T, K, Kneg));
@test_throws(AssertionError, BarrierOptionDownIn(T, Kneg, D));
@test_throws(AssertionError, BarrierOptionDownIn(T, K, Kneg));
@test_throws(AssertionError, BarrierOptionUpIn(T, Kneg, D));
@test_throws(AssertionError, BarrierOptionUpIn(T, K, Kneg));
@test_throws(AssertionError, BarrierOptionUpOut(T, Kneg, D));
@test_throws(AssertionError, BarrierOptionUpOut(T, K, Kneg));
@test_throws(AssertionError, DoubleBarrierOption(T, K, K, Kneg));
@test_throws(AssertionError, DoubleBarrierOption(T, Kneg, K, K));
@test_throws(AssertionError, DoubleBarrierOption(T, K, Kneg, K));

@test FinancialMonteCarlo.maturity(Spot()) == 0.0
@test 2.0 * EuropeanOption(T, K) == EuropeanOption(T, K) + 1.0 * EuropeanOption(T, K)
@test 2.0 * EuropeanOption(T, K) == EuropeanOption(T, K) + EuropeanOption(T, K)

########

import FinancialMonteCarlo.AbstractPayoff

struct tmptype2 <: FinancialMonteCarlo.AbstractPayoff{Float64}
    T::Float64
end

@test_throws(ErrorException, payoff([1.0 1.0; 1.0 1.0], tmptype2(1.0), rfCurve, mc1));
