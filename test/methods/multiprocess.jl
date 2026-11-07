include("../testutils.jl")
using Distributed, ForwardDiff, Random, Statistics, VectorizedRNG

# On real workers, what the serialization of the batches and of their results needs is tested too.
# A batch is a serial run with its own seed, so the results are those of the serial runs.
worker_ids = addprocs(2; exeflags = "--project=$(Base.active_project())")
@everywhere using FinancialMonteCarlo, ForwardDiff
# Differentiated functions must exist on the workers. @everywhere defines them in Main.
@everywhere begin
    mp_config() = MonteCarloConfiguration(4000, 5; mode = FinancialMonteCarlo.MultiProcess(4), seed = [11, 12, 13, 14])
    mp_model(s) = BlackScholesProcess(0.2, Underlying(s, 0.01))
    mp_eu(k) = EuropeanOption(1.0, k)
    mp_by_vol(s) = pricer(BlackScholesProcess(s, Underlying(100.0, 0.01)), ZeroRate(0.02), mp_config(), mp_eu(100.0))
    mp_by_weight(w) = pricer(mp_model(100.0), ZeroRate(0.02), mp_config(), w * mp_eu(100.0) + Forward(1.0))
    mp_by_strike(k) = pricer(mp_model(100.0), ZeroRate(0.02), mp_config(), 2.0 * mp_eu(k) + Forward(1.0))
    mp_by_spot(s) = pricer("A" → mp_model(s), ZeroRate(0.02), mp_config(), "A" → mp_eu(100.0))
end

@testset "MultiProcess" begin
    FMC = FinancialMonteCarlo
    @test nworkers() == 2
    model = BlackScholesProcess(sigma, Underlying(S0, d))
    eu, fw = EuropeanOption(T, K), Forward(T)
    Nstep_ = 5
    seeds = [11, 12, 13, 14]
    process(Nsim_, sampling = FMC.StandardMC(); kwargs...) = MonteCarloConfiguration(Nsim_, Nstep_; sampling = sampling, mode = FMC.MultiProcess(4), seed = seeds, kwargs...)
    mc_mp = process(4000)
    batches = [MonteCarloConfiguration(1000, Nstep_; seed = seed) for seed in seeds]
    on_serial(f) = mean(f(batch) for batch in batches)

    @testset "prices are the average of the serial batches" begin
        smile = FMC.EuropeanOptionSmile(T, [90.0, 110.0], [true, false])
        for payoffs in (eu, [eu, fw], smile, 2.0 * eu + fw)
            @test pricer(model, rfCurve, mc_mp, payoffs) ≈ on_serial(mc -> pricer(model, rfCurve, mc, payoffs))
        end
        @test pricer(model, rfCurve, mc_mp, [eu, fw]) isa Vector{Float64}
        @test pricer(model, rfCurve, mc_mp, smile) isa Vector{Float64}
        @test pricer("A" → model, rfCurve, mc_mp, "A" → eu) ≈ on_serial(mc -> pricer("A" → model, rfCurve, mc, "A" → eu))
        model_nd = GaussianCopulaNVariateProcess(model, model, 0.4)
        eu_nd = EuropeanOptionND(T, K)
        @test pricer(model_nd, rfCurve, mc_mp, eu_nd) ≈ on_serial(mc -> pricer(model_nd, rfCurve, mc, eu_nd))
    end

    @testset "reproducible with $name" for (name, rng) in (("MersenneTwister", MersenneTwister()), ("Xoshiro", Xoshiro()), ("default_rng", Random.default_rng()))
        mc_rng = process(1000; rng = rng)
        prices = [pricer(model, rfCurve, mc_rng, eu) for _ = 1:5]
        @test allequal(prices)
    end

    @testset "metrics" begin
        foreach(FMC.set_seed!, batches)
        S = reduce(vcat, [simulate(model, rfCurve, batch, T) for batch in batches])
        FMC.set_seed!(mc_mp)
        @test simulate(model, rfCurve, mc_mp, T) == S
        @test FMC.distribution(model, rfCurve, mc_mp, T) == S[:, end]
        Payoff = payoff(S, eu, rfCurve, mc_mp)
        @test variance(model, rfCurve, mc_mp, eu) ≈ var(Payoff)
        @test variance(model, rfCurve, mc_mp, [eu, fw])[1] ≈ var(Payoff)
        @test length(variance(model, rfCurve, mc_mp, FMC.EuropeanOptionSmile(T, [90.0, 110.0], [true, false]))) == 2
        american = AmericanOption(T, K, false)
        for payoff_ in (eu, american)
            lower, upper = confinter(model, rfCurve, mc_mp, payoff_)
            @test (lower + upper) / 2 ≈ pricer(model, rfCurve, mc_mp, payoff_)
        end
        @test length(confinter(model, rfCurve, mc_mp, [eu, fw])) == 2
        @test delta(model, rfCurve, mc_mp, eu, 1e-5) ≈ on_serial(mc -> delta(model, rfCurve, mc, eu, 1e-5))
        @test rho(model, rfCurve, mc_mp, eu, 1e-5) ≈ on_serial(mc -> rho(model, rfCurve, mc, eu, 1e-5))
        @test delta(model, rfCurve, mc_mp, 2.0 * eu + fw, 1e-5) ≈ on_serial(mc -> delta(model, rfCurve, mc, 2.0 * eu + fw, 1e-5))
    end

    @testset "paths and seeds" begin
        @test FMC.batch_sizes(process(10)) == [3, 3, 2, 2]
        @test FMC.batch_sizes(process(10, FMC.AntitheticMC())) == [4, 2, 2, 2]
        @test size(simulate(model, rfCurve, process(10, FMC.AntitheticMC()), T), 1) == 10
        @test_throws ArgumentError pricer(model, rfCurve, process(3), eu)
        @test [state.seed for state in MonteCarloConfiguration(10, 1; mode = FMC.MultiProcess(2), seed = 5).rng] == [6, 7]
        @test FMC.MultiProcess().nbatches == nworkers()
        @test_throws AssertionError MonteCarloConfiguration(10, 1; mode = FMC.MultiProcess(2), seed = [1])
        @test_throws AssertionError FMC.MultiProcess(0)
    end

    @testset "control variates" begin
        pilot = MonteCarloConfiguration(1000, Nstep_; seed = 1)
        for (variate, payoff_) in ((fw, eu), (eu, fw))
            mc_cv = process(40_000, FMC.ControlVariates(variate, pilot))
            price = pricer(model, rfCurve, mc_cv, payoff_)
            @test price == pricer(model, rfCurve, mc_cv, payoff_)
            @test abs(price - pricer(model, rfCurve, process(40_000), payoff_)) < 1.0
        end
    end

    @testset "not supported" begin
        @test_throws ArgumentError process(1024, FMC.SobolMC())
        @test_throws ArgumentError process(4000; reseed = false)
        @test_throws ArgumentError MonteCarloConfiguration(10, 1; mode = FMC.MultiProcess(2), rng = local_rng())
        @test_throws MethodError FMC.MultiProcess(FMC.CudaMode())
    end

    @testset "forward automatic differentiation" begin
        @test ForwardDiff.derivative(Main.mp_by_vol, sigma) ≈ on_serial(mc -> ForwardDiff.derivative(s -> pricer(BlackScholesProcess(s, Underlying(S0, d)), rfCurve, mc, eu), sigma))
        @test ForwardDiff.derivative(Main.mp_by_weight, 2.0) ≈ on_serial(mc -> ForwardDiff.derivative(w -> pricer(model, rfCurve, mc, w * eu + fw), 2.0))
        @test ForwardDiff.derivative(Main.mp_by_strike, K) ≈ on_serial(mc -> ForwardDiff.derivative(k -> pricer(model, rfCurve, mc, 2.0 * EuropeanOption(T, k) + fw), K))
        @test ForwardDiff.derivative(Main.mp_by_spot, S0) ≈ on_serial(mc -> ForwardDiff.derivative(s -> pricer("A" → BlackScholesProcess(sigma, Underlying(s, d)), rfCurve, mc, "A" → eu), S0))
        # the workers do not have this function: the error must not become a wrong derivative
        @test_throws Exception ForwardDiff.derivative(s -> pricer(BlackScholesProcess(s, Underlying(S0, d)), rfCurve, mc_mp, eu), sigma)
    end

    @testset "the number of workers does not matter" begin
        price = pricer(model, rfCurve, mc_mp, eu)
        rmprocs(worker_ids[2])
        @test pricer(model, rfCurve, mc_mp, eu) == price
    end
end

rmprocs(workers())
