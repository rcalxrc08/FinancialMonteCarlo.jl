include("../testutils.jl")
using Distributed, ForwardDiff, Random, Statistics, VectorizedRNG

@testset "Parallel modes and RNGs" begin
    model = BlackScholesProcess(sigma, Underlying(S0, d))
    eu = EuropeanOption(T, K)

    @testset "$name" for (name, mode) in [
        ("MultiThreading", FinancialMonteCarlo.MultiThreading()),
        ("MultiProcess", FinancialMonteCarlo.MultiProcess()),
        ("VectorizedRNG", FinancialMonteCarlo.SerialMode(10, local_rng())),
    ]
        mc_par = MonteCarloConfiguration(Nsim, Nstep, mode)
        @test abs(pricer(model, rfCurve, mc_par, eu) - BS_CALL_REF.eu) < 1.0
        prices = pricer(model, rfCurve, mc_par, [eu, eu])
        @test prices[1] ≈ prices[2]
    end
end

# A batch is a serial run with its own seed, so on any number of threads the results are those of the serial runs.
@testset "MultiThreading" begin
    FMC = FinancialMonteCarlo
    model = BlackScholesProcess(sigma, Underlying(S0, d))
    eu, fw = EuropeanOption(T, K), Forward(T)
    Nstep_ = 5
    seeds = [11, 12, 13, 14]
    threaded(Nsim_, method = FMC.StandardMC(), args...) = MonteCarloConfiguration(Nsim_, Nstep_, method, FMC.MultiThreading(FMC.SerialMode(), 4, seeds), args...)
    mc_mt = threaded(4000)
    batches = [MonteCarloConfiguration(1000, Nstep_, FMC.StandardMC(), FMC.SerialMode(seed)) for seed in seeds]
    on_serial(f) = mean(f(batch) for batch in batches)
    function serial_paths()
        foreach(FMC.set_seed!, batches)
        return reduce(vcat, [simulate(model, rfCurve, batch, T) for batch in batches])
    end

    @testset "prices are the average of the serial batches" begin
        smile = FMC.EuropeanOptionSmile(T, [90.0, 110.0], [true, false])
        for payoffs in (eu, [eu, fw], smile, 2.0 * eu + fw)
            @test pricer(model, rfCurve, mc_mt, payoffs) ≈ on_serial(mc -> pricer(model, rfCurve, mc, payoffs))
        end
        @test pricer(model, rfCurve, mc_mt, [eu, fw]) isa Vector{Float64}
        @test pricer(model, rfCurve, mc_mt, smile) isa Vector{Float64}
        @test pricer("A" → model, rfCurve, mc_mt, "A" → eu) ≈ on_serial(mc -> pricer("A" → model, rfCurve, mc, "A" → eu))
    end

    @testset "multivariate" begin
        model_nd = GaussianCopulaNVariateProcess(model, model, 0.4)
        eu_nd = EuropeanOptionND(T, K)
        @test pricer(model_nd, rfCurve, mc_mt, eu_nd) ≈ on_serial(mc -> pricer(model_nd, rfCurve, mc, eu_nd))
        @test pricer("A_B" → model_nd, rfCurve, mc_mt, "A_B" → eu_nd) ≈ on_serial(mc -> pricer("A_B" → model_nd, rfCurve, mc, "A_B" → eu_nd))
        FMC.set_seed!(mc_mt)
        paths_mt = simulate(model_nd, rfCurve, mc_mt, T)
        foreach(FMC.set_seed!, batches)
        paths = [simulate(model_nd, rfCurve, batch, T) for batch in batches]
        @test length(paths_mt) == 2
        @test all(paths_mt[i] == reduce(vcat, [batch[i] for batch in paths]) for i = 1:2)
        @test variance(model_nd, rfCurve, mc_mt, eu_nd) ≈ var(payoff(paths_mt, eu_nd, rfCurve, mc_mt))
        positions = Array{FMC.Position}(undef, 3)
        positions[1], positions[3] = 1.0 * EuropeanOption(T, K), 1.0 * eu_nd
        @test pricer(model_nd, rfCurve, mc_mt, positions) ≈ on_serial(mc -> pricer(model_nd, rfCurve, mc, positions))
    end

    @testset "the generators are not shared" begin
        rngs = [mode.rng for mode in mc_mt.parallelMode.sub_modes]
        @test allunique(objectid.(rngs))
        @test allunique(pricer(model, rfCurve, MonteCarloConfiguration(1000, Nstep_, FMC.StandardMC(), FMC.SerialMode(seed)), eu) for seed in seeds)
    end

    @testset "reproducible with $name" for (name, rng) in (("MersenneTwister", MersenneTwister()), ("Xoshiro", Xoshiro()), ("default_rng", Random.default_rng()))
        mc_rng = MonteCarloConfiguration(1000, Nstep_, FMC.MultiThreading(FMC.SerialMode(0, rng), 4, seeds))
        prices = [pricer(model, rfCurve, mc_rng, eu) for _ = 1:50]
        @test allequal(prices)
    end

    @testset "metrics on the paths of all the batches" begin
        S = serial_paths()
        @test simulate(model, rfCurve, (FMC.set_seed!(mc_mt); mc_mt), T) == S
        @test FMC.distribution(model, rfCurve, mc_mt, T) == S[:, end]
        Payoff = payoff(S, eu, rfCurve, mc_mt)
        @test variance(model, rfCurve, mc_mt, eu) ≈ var(Payoff)
        @test variance(model, rfCurve, mc_mt, [eu, fw])[1] ≈ var(Payoff)
        for payoff_ in (eu, AmericanOption(T, K, false))
            lower, upper = confinter(model, rfCurve, mc_mt, payoff_)
            @test (lower + upper) / 2 ≈ pricer(model, rfCurve, mc_mt, payoff_)
        end
        @test length(confinter(model, rfCurve, mc_mt, [eu, fw])) == 2
        @test delta(model, rfCurve, mc_mt, eu, 1e-5) ≈ on_serial(mc -> delta(model, rfCurve, mc, eu, 1e-5))
        @test rho(model, rfCurve, mc_mt, eu, 1e-5) ≈ on_serial(mc -> rho(model, rfCurve, mc, eu, 1e-5))
        @test delta(model, rfCurve, mc_mt, 2.0 * eu + fw, 1e-5) ≈ on_serial(mc -> delta(model, rfCurve, mc, 2.0 * eu + fw, 1e-5))
        @test size(simulate(BrownianMotion(sigma, 0.01), mc_mt, T)) == (4000, Nstep_ + 1)
    end

    @testset "paths" begin
        @test FMC.batch_sizes(threaded(10)) == [3, 3, 2, 2]
        @test FMC.batch_sizes(threaded(10, FMC.AntitheticMC())) == [4, 2, 2, 2]
        @test size(simulate(model, rfCurve, threaded(10), T), 1) == 10
        @test size(simulate(model, rfCurve, threaded(10, FMC.AntitheticMC()), T), 1) == 10
        @test abs(pricer(model, rfCurve, threaded(40_002), eu) - BS_CALL_REF.eu) < 1.0
        @test_throws ArgumentError pricer(model, rfCurve, threaded(3), eu)
        @test_throws ArgumentError pricer(model, rfCurve, threaded(6, FMC.AntitheticMC()), eu)
    end

    @testset "seeds" begin
        @test [mode.seed for mode in FMC.MultiThreading(FMC.SerialMode(5), 2).sub_modes] == [6, 7]
        @test length(FMC.MultiThreading(FMC.SerialMode(), Threads.nthreads() + 1).sub_modes) == Threads.nthreads() + 1
        @test_throws AssertionError FMC.MultiThreading(FMC.SerialMode(), 2, [1])
        @test_throws AssertionError FMC.MultiThreading(FMC.SerialMode(), 0)
    end

    @testset "init_rng = false continues the streams of the batches" begin
        mc_stream = threaded(4000, FMC.StandardMC(), false)
        first_price = pricer(model, rfCurve, mc_stream, eu)
        @test first_price ≈ pricer(model, rfCurve, mc_mt, eu)
        @test pricer(model, rfCurve, mc_stream, eu) != first_price
    end

    @testset "control variates" begin
        pilot = MonteCarloConfiguration(1000, Nstep_, FMC.StandardMC(), FMC.SerialMode(1))
        for (variate, payoff_) in ((fw, eu), (eu, fw))
            mc_cv = MonteCarloConfiguration(40_000, Nstep_, FMC.ControlVariates(variate, pilot), FMC.MultiThreading(FMC.SerialMode(), 4, seeds))
            price = pricer(model, rfCurve, mc_cv, payoff_)
            @test price == pricer(model, rfCurve, mc_cv, payoff_)
            @test abs(price - pricer(model, rfCurve, threaded(40_000), payoff_)) < 1.0
        end
        cv_method = FMC.ControlVariates(fw, pilot)
        price = pricer(model, rfCurve, threaded(4000, cv_method), eu)
        mc_stream = threaded(4000, cv_method, false)
        @test pricer(model, rfCurve, mc_stream, eu) ≈ price
        @test pricer(model, rfCurve, mc_stream, eu) != price
    end

    @testset "not supported" begin
        @test_throws ArgumentError pricer(model, rfCurve, threaded(1024, FMC.SobolMode()), eu)
        @test_throws MethodError FMC.MultiThreading(FMC.CudaMode())
        @test_throws MethodError FMC.MultiThreading(FMC.MultiThreading())
        @test_throws ArgumentError FMC.MultiThreading(FMC.SerialMode(0, local_rng()), 2)
    end

    @testset "forward automatic differentiation" begin
        by_vol(s, mc) = pricer(BlackScholesProcess(s, Underlying(S0, d)), rfCurve, mc, eu)
        @test ForwardDiff.derivative(s -> by_vol(s, mc_mt), sigma) ≈ on_serial(mc -> ForwardDiff.derivative(s -> by_vol(s, mc), sigma))
        by_weight(w, mc) = pricer(model, rfCurve, mc, w * eu + fw)
        @test ForwardDiff.derivative(w -> by_weight(w, mc_mt), 2.0) ≈ on_serial(mc -> ForwardDiff.derivative(w -> by_weight(w, mc), 2.0))
        by_strike(k, mc) = pricer(model, rfCurve, mc, 2.0 * EuropeanOption(T, k) + fw)
        @test ForwardDiff.derivative(k -> by_strike(k, mc_mt), K) ≈ on_serial(mc -> ForwardDiff.derivative(k -> by_strike(k, mc), K))
        by_spot(s, mc) = pricer("A" → BlackScholesProcess(sigma, Underlying(s, d)), rfCurve, mc, "A" → eu)
        @test ForwardDiff.derivative(s -> by_spot(s, mc_mt), S0) ≈ on_serial(mc -> ForwardDiff.derivative(s -> by_spot(s, mc), S0))
    end
end

# CudaMode exists without CUDA, so that the package extensions of other packages can dispatch on it.
@testset "CudaMode without CUDA" begin
    mode = FinancialMonteCarlo.CudaMode(1, 2)
    @test mode isa FinancialMonteCarlo.AbstractCudaMode
    @test (mode.seed_cpu, mode.seed_gpu) == (1, 2)
end
