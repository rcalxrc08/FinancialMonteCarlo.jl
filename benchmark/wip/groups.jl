#Number Type (just sigma tested) ---> Model type ----> Mode ---> zero rate type

using BenchmarkTools, DualNumbers, FinancialMonteCarlo, VectorizedRNG

const sigma_dual = dual(0.2, 1.0);
const sigma_no_dual = 0.2;

suite_num = BenchmarkGroup()

const und = Underlying(100.0, 0.01);
const p = 0.3;
const lam = 5.0;
const lamp = 30.0;
const lamm = 20.0;
mu1 = 0.03;
sigma1 = 0.02;

sigma_zero = 0.2;
kappa = 0.01;
theta = 0.03;
lambda = 0.01;
rho = 0.0;

theta1 = 0.01;
k1 = 0.03;
sigma1 = 0.02;

T = 1.0;
K = 100.0;

const Nsim = 10000;
const Nstep = 30;
const std_mc = MonteCarloConfiguration(Nsim, Nstep);
const anti_mc = MonteCarloConfiguration(Nsim, Nstep, FinancialMonteCarlo.AntitheticMC());
const sobol_mc = MonteCarloConfiguration(Nsim, Nstep, FinancialMonteCarlo.SobolMode());
const cv_mc = MonteCarloConfiguration(Nsim, Nstep, FinancialMonteCarlo.ControlVariates(Forward(T), MonteCarloConfiguration(100, 3)));
const vect_mc = MonteCarloConfiguration(Nsim, Nstep, FinancialMonteCarlo.SerialMode(), 10, local_rng());

const r = 0.02;
const r_prev = 0.019999;
lam_vec = Float64[0.999999999];
const zr_scalar = ZeroRate(r);
const zr_imp = FinancialMonteCarlo.ImpliedZeroRate([r_prev, r], 2.0);

const opt_ = EuropeanOption(T, K);

bs(sigma) = BlackScholesProcess(sigma, und)
kou(sigma) = KouProcess(sigma, lam, p, lamp, lamm, und);
hest(sigma) = HestonProcess(sigma, sigma_zero, lambda, kappa, rho, theta, und);
vg(sigma) = VarianceGammaProcess(sigma, theta1, k1, und)
merton(sigma) = MertonProcess(sigma, lambda, mu1, sigma1, und)
nig(sigma) = NormalInverseGaussianProcess(sigma, theta1, k1, und);
shift_mixture(sigma) = ShiftedLogNormalMixture([sigma, 0.2], lam_vec, 0.0, und);

for sigma in Number[sigma_no_dual, sigma_dual]
    for model_ in [bs(sigma), kou(sigma), hest(sigma), vg(sigma), merton(sigma), nig(sigma), shift_mixture(sigma)]
        for mode_ in [std_mc, anti_mc, sobol_mc, cv_mc, vect_mc]
            for zero_ in [zr_scalar, zr_imp]
                suite_num[string(typeof(sigma)), string(typeof(model_).name), string(typeof(mode_.monteCarloMethod)), string(typeof(mode_.rng)), string(typeof(zero_))] = @benchmarkable pricer($model_, $zero_, $mode_, $opt_)
            end
        end
    end
end
#loadparams!(suite, BenchmarkTools.load("params.json")[1], :evals, :samples);
results = run(suite_num, verbose = true)