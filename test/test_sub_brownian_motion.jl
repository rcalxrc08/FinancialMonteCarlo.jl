using FinancialMonteCarlo, Distributions, Test;
@show "Test Parameters"
S0 = 100.0;
r = 0.02;
d = 0.01;
T = 1.0;
Nsim = 100;
Nstep = 30;

sigma = 0.2;

McConfig = MonteCarloConfiguration(Nsim, Nstep);
mc1 = MonteCarloConfiguration(Nsim, Nstep, FinancialMonteCarlo.AntitheticMC());
rfCurve = ZeroRate(r);
rfCurve2 = ZeroRate([0.00, 0.02], T);

@show "Test Subordinated Brownian Motion Parameters"
drift = 0.0;
@test_throws(ErrorException, SubordinatedBrownianMotion(-sigma, drift, InverseGaussian(1.0, 1.0)))
@test_throws(ErrorException, SubordinatedBrownianMotion(-sigma, rfCurve2.r, InverseGaussian(1.0, 1.0)))
