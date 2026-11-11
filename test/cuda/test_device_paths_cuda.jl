# What runs on the device instead of the host in CudaMode (#29): the jumps added by one kernel, the gamma and inverse
# Gaussian subordinators drawn on the device, barrier and Asian payoffs as row reductions. Needs a GPU, not run on CI.
using Test, FinancialMonteCarlo, CUDA, Distributions, Statistics
const FMC = FinancialMonteCarlo
const Ext = Base.get_extension(FMC, :FinancialMonteCarloCUDAExt)

S0, K, r, T, d, sigma, D = 100.0, 100.0, 0.02, 1.0, 0.01, 0.2, 90.0
rfCurve = ZeroRate(r)
und = Underlying(S0, d)
Nsim, Nstep = 100_000, 30
gpu(sampling = FMC.StandardMC()) = MonteCarloConfiguration(Nsim, Nstep; sampling = sampling, mode = FMC.CudaMode(), seed = 1)
cpu(sampling = FMC.StandardMC()) = MonteCarloConfiguration(Nsim, Nstep; sampling = sampling, seed = 3)

@testset "device subordinators" begin
    mc = gpu()
    for law in (Gamma(0.02 / 0.03, 0.03), Gamma(2.5, 0.4), InverseGaussian(0.02, 0.02^2 / 0.03), InverseGaussian(1.0, 3.0))
        x = CuArray{Float64}(undef, 1_000_000)
        FMC.set_seed!(mc)
        h = Array(Ext.device_draws!(x, law, mc))
        FMC.set_seed!(mc)
        @test Array(Ext.device_draws!(similar(x), law, mc)) == h # reproducible
        @test all(>=(0), h)
        @test abs(mean(h) - mean(law)) < 4 * std(law) / sqrt(length(h))
        @test abs(var(h) / var(law) - 1) < 0.02
    end
    # Float32
    x32 = Ext.device_draws!(CuArray{Float32}(undef, 1000), Gamma(0.5f0, 2.0f0), gpu())
    @test eltype(x32) == Float32 && all(isfinite, Array(x32))
end

@testset "subordinated models: GPU as CPU" begin
    for model in (VarianceGammaProcess(sigma, 0.01, 0.03, und), NormalInverseGaussianProcess(sigma, 0.01, 0.03, und)), sampling in (FMC.StandardMC(), FMC.AntitheticMC())
        @test abs(pricer(model, rfCurve, gpu(sampling), EuropeanOption(T, K)) - pricer(model, rfCurve, cpu(sampling), EuropeanOption(T, K))) < 0.15
    end
end

@testset "jumps with one kernel: GPU as CPU" begin
    for model in (KouProcess(sigma, 3.0, 0.5, 25.0, 30.0, und), MertonProcess(sigma, 3.0, 0.03, 0.02, und)), sampling in (FMC.StandardMC(), FMC.AntitheticMC())
        # the same host draws: the jumps of the two simulations are the same
        X = CUDA.zeros(Float64, 1000, Nstep + 1)
        mc_small = MonteCarloConfiguration(1000, Nstep; sampling = sampling, mode = FMC.CudaMode(), seed = 5)
        FMC.set_seed!(mc_small)
        FMC.add_jump_matrix!(X, model, mc_small, T)
        Y = zeros(1000, Nstep + 1)
        FMC.set_seed!(mc_small)
        FMC.add_jumps!(Y, model, mc_small, T)
        @test Array(X) ≈ Y
        @test abs(pricer(model, rfCurve, gpu(sampling), EuropeanOption(T, K)) - pricer(model, rfCurve, cpu(sampling), EuropeanOption(T, K))) < 0.15
    end
end

@testset "barrier and Asian payoffs on the device" begin
    model = BlackScholesProcess(sigma, und)
    mc = gpu()
    S = FMC.simulate(model, rfCurve, mc, T)
    for p in (BarrierOptionDownOut(T, K, D), BarrierOptionDownIn(T, K, D), BarrierOptionUpOut(T, K, 120.0), BarrierOptionUpIn(T, K, 120.0), DoubleBarrierOption(T, K, D, 120.0), AsianFixedStrikeOption(T, K), AsianFloatingStrikeOption(T))
        v = FMC.payoff(S, p, rfCurve, mc)
        @test v isa CuVector{Float64}
        @test Array(v) ≈ FMC.payoff(Array(S), p, rfCurve, mc)
        @test abs(pricer(model, rfCurve, mc, p) - pricer(model, rfCurve, cpu(), p)) < 0.15
    end
end

@testset "other subordinator laws: on the host, copied" begin
    # a law with no device sampler: drawn by Distributions on the host
    model = SubordinatedBrownianMotion(sigma, 0.01, 0.1, Exponential(T / Nstep))
    S = simulate(model, gpu(), T)
    @test S isa CuMatrix{Float64}
    @test abs(mean(Array(S)[:, end]) - mean(simulate(model, cpu(), T)[:, end])) < 0.01
end
