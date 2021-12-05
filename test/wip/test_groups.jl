#Number Type (just sigma tested) ---> Model type ----> Mode ---> zero rate type

using Test, DualNumbers, FinancialMonteCarlo, VectorizedRNG, JLD2
rebase = false;
sigma_dual = dual(0.2, 1.0);
sigma_no_dual = 0.2;

suite_num = Dict{String, Number}()

und = Underlying(100.0, 0.01);
p = 0.3;
lam = 5.0;
lamp = 30.0;
lamm = 20.0;
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

Nsim = 10_000;
Nstep = 30;
std_mc = MonteCarloConfiguration(Nsim, Nstep);
anti_mc = MonteCarloConfiguration(Nsim, Nstep, FinancialMonteCarlo.AntitheticMC());
sobol_mc = MonteCarloConfiguration(Nsim, Nstep, FinancialMonteCarlo.SobolMode());
cv_mc = MonteCarloConfiguration(Nsim, Nstep, FinancialMonteCarlo.ControlVariates(Forward(T), MonteCarloConfiguration(100, 3)));
vect_mc = MonteCarloConfiguration(Nsim, Nstep, FinancialMonteCarlo.SerialMode(), 10, local_rng());

r = 0.02;
r_prev = 0.019999;
D = 60;
lam_vec = Float64[0.99];
zr_scalar = ZeroRate(r);
zr_imp = FinancialMonteCarlo.ImpliedZeroRate([r_prev, r], 2.0);

eur_opt_ = EuropeanOption(T, K);
EUBin = BinaryEuropeanOption(T, K)
AMData = AmericanOption(T, K)
AmBin = BinaryEuropeanOption(T, K)
BarrierData = BarrierOptionDownOut(T, K, D)
AsianFloatingStrikeData = AsianFloatingStrikeOption(T)
AsianFixedStrikeData = AsianFixedStrikeOption(T, K)
BinAMData = BinaryAmericanOption(T, K)
BarrierDataDI = BarrierOptionDownIn(T, K, D)
BarrierDataUI = BarrierOptionUpIn(T, K, D)
BarrierDataUO = BarrierOptionUpOut(T, K, D)
BermData = BermudanOption(collect(0.2:0.1:T), K);
doubleBarrierOptionDownOut = DoubleBarrierOption(T, K, K / 10.0, 1.2 * K)

bs(sigma) = BlackScholesProcess(sigma, und)
kou(sigma) = KouProcess(sigma, lam, p, lamp, lamm, und);
hest(sigma) = HestonProcess(sigma, sigma_zero, lambda, kappa, rho, theta, und);
vg(sigma) = VarianceGammaProcess(sigma, theta1, k1, und)
merton(sigma) = MertonProcess(sigma, lambda, mu1, sigma1, und)
nig(sigma) = NormalInverseGaussianProcess(sigma, theta1, k1, und);
log_mixture(sigma) = LogNormalMixture([sigma, 0.2], lam_vec, und);
shift_mixture(sigma) = ShiftedLogNormalMixture([sigma, 0.2], lam_vec, 0.0, und);
create_str_lam(x...) = join([string(typeof(y)) for y in x], "_");

for sigma in Number[sigma_no_dual, sigma_dual]
    for model_ in [bs(sigma), kou(sigma), hest(sigma), vg(sigma), merton(sigma), nig(sigma), log_mixture(sigma), shift_mixture(sigma)]
        for mode_ in [std_mc, anti_mc, sobol_mc, cv_mc, vect_mc]
            for zero_ in [zr_scalar, zr_imp]
                for opt in [eur_opt_, EUBin, AMData, AmBin, BarrierData, AsianFloatingStrikeData, AsianFixedStrikeData, BinAMData, BarrierDataDI, BarrierDataUO, BermData, doubleBarrierOptionDownOut]
                    key_tmp1 = create_str_lam(sigma, model_, mode_.monteCarloMethod, mode_.rng, zero_, opt)
                    result = pricer(model_, zero_, mode_, opt)
                    @assert !isnan(result) "Result for provided simulation is nan $(result)  $(key_tmp1)\n"
                    suite_num[key_tmp1] = result
                end
            end
        end
    end
end

if rebase
    save("test_results.jld2", suite_num)
end
suite_num_old_results = load("test_results.jld2")
toll = 1.0

for key_tmp in unique([collect(keys(suite_num_old_results)); collect(keys(suite_num))])
    @test abs(suite_num[key_tmp] - suite_num_old_results[key_tmp]) < toll
end