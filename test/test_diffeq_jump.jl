using Test, FinancialMonteCarlo, DifferentialEquations, Random;
@show "Differential Equation Junctor"
Nsim = 10000;
Nstep = 30;
toll = 0.8;
mc = MonteCarloConfiguration(Nsim, Nstep);
mc1 = MonteCarloConfiguration(Nsim, Nstep, FinancialMonteCarlo.AntitheticMC());
S0 = 100.0
K = 100.0;
D = 90.0;
r = 0.02
sigma = 0.2
T = 1.0;
d = 0.01;
u0 = 0.0;

p1 = 0.3;
lam = 5.0;
lamp = 30.0;
lamm = 20.0;

#Drift
dr_(u, p, t) = (r - d - sigma^2 / 2.0 - lam * (p1 / (lamp - 1) - (1 - p1) / (lamm + 1)))
#Diffusion
g_1(u, p, t) = sigma
#Time Window
tspan = (0.0, T)
rate(u, p, t) = lam * T
affect!(integrator) = (integrator.u = integrator.u + ((rand(mc.parallelMode.rng) < p1) ? randexp(mc.parallelMode.rng) / lamp : -randexp(mc.parallelMode.rng) / lamm))
jump = ConstantRateJump(rate, affect!)

#Definition of the SDE
prob = SDEProblem(dr_, g_1, u0, tspan)
jump_prob = JumpProblem(prob, Direct(), jump, rng = mc.parallelMode.rng)
monte_prob = EnsembleProblem(jump_prob)
rfCurve = ZeroRate(r);
func(x) = S0 * exp(x);
model = FinancialMonteCarlo.MonteCarloDiffEqModel(monte_prob, func, Underlying(S0, d))

FwdData = Forward(T)
EUData = EuropeanOption(T, K)
AMData = AmericanOption(T, K)
BarrierData = BarrierOptionDownOut(T, K, D)
AsianFloatingStrikeData = AsianFloatingStrikeOption(T)
AsianFixedStrikeData = AsianFixedStrikeOption(T, K)

@show FwdPrice = pricer(model, rfCurve, mc, FwdData);
@show EuPrice = pricer(model, rfCurve, mc, EUData);
@show AmPrice = pricer(model, rfCurve, mc, AMData);
@show BarrierPrice = pricer(model, rfCurve, mc, BarrierData);
@show AsianPrice1 = pricer(model, rfCurve, mc, AsianFloatingStrikeData);
@show AsianPrice2 = pricer(model, rfCurve, mc, AsianFixedStrikeData);

@test abs(FwdPrice - 98.8436678850961) < toll

Tneg = -T;
tspanNeg = (0.0, Tneg)
probNeg = SDEProblem(dr_, g_1, u0, tspanNeg, rng = mc.parallelMode.rng)
monte_probNeg = EnsembleProblem(probNeg)
