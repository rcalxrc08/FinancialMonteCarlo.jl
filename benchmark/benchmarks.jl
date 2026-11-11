# The benchmark suite: `SUITE`, a BenchmarkTools group of groups, each benchmark one call of a metric.
# Run it with benchmark/run.jl, which prints the table of the results; PkgBenchmark and AirspeedVelocity also read `SUITE`.
using BenchmarkTools, FinancialMonteCarlo, ForwardDiff
const FMC = FinancialMonteCarlo

const SUITE = BenchmarkGroup()

# Market data and simulation settings
const S0, K, r, T, d, D, σ = 100.0, 100.0, 0.02, 1.0, 0.01, 90.0, 0.2
const Nsim, Nstep = 10_000, 50
const underlying = Underlying(S0, d)
const zero_rate = ZeroRate(r)
const mc = MonteCarloConfiguration(Nsim, Nstep)
const european = EuropeanOption(T, K)

const models = (BlackScholes = BlackScholesProcess(σ, underlying), Heston = HestonProcess(σ, 0.2, 1.5, 0.04, -0.7, 0.04, underlying), Kou = KouProcess(σ, 3.0, 0.5, 25.0, 30.0, underlying), Merton = MertonProcess(σ, 3.0, -0.05, 0.1, underlying), VarianceGamma = VarianceGammaProcess(σ, -0.1, 0.2, underlying), NIG = NormalInverseGaussianProcess(σ, -0.1, 0.2, underlying), LogNormalMixture = LogNormalMixture([0.2, 0.3], [0.4], underlying))

# Every model on a European option
SUITE["models"] = BenchmarkGroup()
for (name, model) in pairs(models)
    SUITE["models"][string(name)] = @benchmarkable pricer($model, $zero_rate, $mc, $european)
end

# The samplings and the execution modes, on the Black-Scholes and Kou models
SUITE["methods"] = BenchmarkGroup()
const methods = (standard = mc, antithetic = MonteCarloConfiguration(Nsim, Nstep; sampling = FMC.AntitheticMC()), sobol = MonteCarloConfiguration(Nsim, Nstep; sampling = FMC.SobolMC()), threads = MonteCarloConfiguration(Nsim, Nstep; mode = FMC.MultiThreading(Threads.nthreads())))
for name in (:BlackScholes, :Kou), (method_name, method) in pairs(methods)
    SUITE["methods"]["$name $method_name"] = @benchmarkable pricer($(models[name]), $zero_rate, $method, $european)
end

# Every family of payoffs on the Black-Scholes model
SUITE["payoffs"] = BenchmarkGroup()
const payoffs = (Forward = Forward(T), European = european, BinaryEuropean = BinaryEuropeanOption(T, K), American = AmericanOption(T, K, false), Bermudan = BermudanOption(collect(0.25:0.25:1.0), K, false), BarrierDownOut = BarrierOptionDownOut(T, K, D), DoubleBarrier = DoubleBarrierOption(T, K, D, 120.0), AsianFixedStrike = AsianFixedStrikeOption(T, K), AsianFloatingStrike = AsianFloatingStrikeOption(T))
for (name, payoff_) in pairs(payoffs)
    SUITE["payoffs"][string(name)] = @benchmarkable pricer($(models.BlackScholes), $zero_rate, $mc, $payoff_)
end
const strikes = collect(80.0:5.0:120.0)
const smile = FMC.EuropeanOptionSmile(T, strikes, fill(true, length(strikes)))
SUITE["payoffs"]["Smile of 9 strikes"] = @benchmarkable pricer($(models.BlackScholes), $zero_rate, $mc, $smile)
const position = EuropeanOption(T, 90.0) - 2 * EuropeanOption(T, K) + EuropeanOption(T, 110.0)
SUITE["payoffs"]["Butterfly position"] = @benchmarkable pricer($(models.BlackScholes), $zero_rate, $mc, $position)

# Rate curves, a basket on a Gaussian copula and a portfolio
SUITE["markets"] = BenchmarkGroup()
const curve_rate = ZeroRateCurve(collect(range(0.01, 0.03; length = 12)), T)
SUITE["markets"]["BlackScholes on a rate curve"] = @benchmarkable pricer($(models.BlackScholes), $curve_rate, $mc, $european)
const copula = GaussianCopulaNVariateProcess(models.BlackScholes, models.Kou, 0.4)
const basket = EuropeanOptionND(T, K)
SUITE["markets"]["Basket on a copula"] = @benchmarkable pricer($copula, $zero_rate, $mc, $basket)
const market = ("A" → models.BlackScholes) + ("B" → models.Heston) + (("C", "D") → GaussianCopulaNVariateProcess(models.BlackScholes, models.Kou, 0.4))
const portfolio = ("A" → european) + ("B" → AmericanOption(T, K, false)) + (("C", "D") → basket)
SUITE["markets"]["Portfolio on four underlyings"] = @benchmarkable pricer($market, $zero_rate, $mc, $portfolio)

# Statistics and sensitivities: bump and reprice, and forward mode AD
SUITE["greeks"] = BenchmarkGroup()
SUITE["greeks"]["variance"] = @benchmarkable variance($(models.BlackScholes), $zero_rate, $mc, $european)
SUITE["greeks"]["confinter"] = @benchmarkable confinter($(models.BlackScholes), $zero_rate, $mc, $european)
SUITE["greeks"]["delta"] = @benchmarkable delta($(models.BlackScholes), $zero_rate, $mc, $european, 1e-4)
SUITE["greeks"]["rho"] = @benchmarkable rho($(models.BlackScholes), $zero_rate, $mc, $european, 1e-4)
bs_price(x) = pricer(BlackScholesProcess(x[1], Underlying(x[2], d)), ZeroRate(x[3]), mc, european)
SUITE["greeks"]["ForwardDiff BlackScholes (σ, S0, r)"] = @benchmarkable ForwardDiff.gradient(bs_price, $([σ, S0, r]))
kou_price(x) = pricer(KouProcess(x[1], x[2], 0.5, 25.0, 30.0, Underlying(S0, d)), zero_rate, mc, european)
SUITE["greeks"]["ForwardDiff Kou (σ, λ)"] = @benchmarkable ForwardDiff.gradient(kou_price, $([σ, 3.0]))
am_price(x) = pricer(BlackScholesProcess(x[1], Underlying(x[2], d)), zero_rate, mc, AmericanOption(T, K, false))
SUITE["greeks"]["ForwardDiff American (σ, S0)"] = @benchmarkable ForwardDiff.gradient(am_price, $([σ, S0]))
