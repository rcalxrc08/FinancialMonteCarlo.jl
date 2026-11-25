using Test, FinancialMonteCarlo, CUDA, ForwardDiff;
const FMC = FinancialMonteCarlo

S0 = 100.0;
r = 0.02;
T = 1.0;
d = 0.01;
sigma = 0.2;

Nsim = 100000;
Nstep = 10;
rfCurve = ZeroRate(r);
model = BlackScholesProcess(sigma);
und = Underlying(S0, d);
K = [80.0, 90.0, 100.0, 110.0, 120.0];
isCall = [true, false, true, false, true];
mc = MonteCarloConfiguration(Nsim, Nstep; mode = FMC.CudaMode(), seed = 1);
mc_cpu = FMC.with_mode(mc, FMC.SerialMode());

@testset "$(nameof(Smile))" for Smile in (FMC.EuropeanOptionSmile, FMC.BinaryEuropeanOptionSmile)
    smile = Smile(T, K, isCall)
    options = FMC.smile_to_options(smile)
    @test which(FMC.payoff, Tuple{CuMatrix{Float64}, typeof(smile), typeof(rfCurve), typeof(mc), Float64}).module === Base.get_extension(FMC, :FinancialMonteCarloCUDAExt)

    @testset "the results stay on the GPU" begin
        @test pricer(model, und, rfCurve, mc, smile) isa CuVector{Float64}
        @test variance(model, und, rfCurve, mc, smile) isa CuVector{Float64}
        @test delta(model, und, rfCurve, mc, smile, 1e-4) isa CuVector{Float64}
        @test rho(model, und, rfCurve, mc, smile, 1e-4) isa CuVector{Float64}
        FMC.set_seed!(mc)
        @test payoff(simulate(model, und, rfCurve, mc, T), smile, rfCurve, mc) isa CuMatrix{Float64}
    end

    @testset "same results of the single options" begin
        @test Array(pricer(model, und, rfCurve, mc, smile)) ≈ pricer(model, und, rfCurve, mc, options)
        @test Array(variance(model, und, rfCurve, mc, smile)) ≈ variance(model, und, rfCurve, mc, options)
        @test all(collect.(confinter(model, und, rfCurve, mc, smile)) .≈ collect.(confinter(model, und, rfCurve, mc, options)))
        FMC.set_seed!(mc)
        S = simulate(model, und, rfCurve, mc, 2 * T)
        @test Array(payoff(S, smile, rfCurve, mc, 2 * T))[:, 2] ≈ Array(payoff(S, options[2], rfCurve, mc, 2 * T))
    end

    @testset "strikes and flags of other array types" begin
        @test pricer(model, und, rfCurve, mc, Smile(T, K, BitVector(isCall))) ≈ pricer(model, und, rfCurve, mc, smile)
        @test pricer(model, und, rfCurve, mc, Smile(T, CuArray(K), CuArray(isCall))) ≈ pricer(model, und, rfCurve, mc, smile)
    end

    # The mode decides where the smile is priced: its strikes move there, once per metric.
    @testset "strikes on the other device" begin
        smile_gpu = Smile(T, CuArray(K), CuArray(isCall))
        @test FMC.to_mode(smile, FMC.CudaMode()).K isa CuVector{Float64}
        @test FMC.to_mode(smile_gpu, FMC.SerialMode()).K isa Vector{Float64}
        @test pricer(model, und, rfCurve, mc_cpu, smile_gpu) isa Vector{Float64}
        @test pricer(model, und, rfCurve, mc_cpu, smile_gpu) ≈ pricer(model, und, rfCurve, mc_cpu, smile)
        @test variance(model, und, rfCurve, mc_cpu, smile_gpu) ≈ variance(model, und, rfCurve, mc_cpu, smile)
    end

    @testset "Float32 stays Float32" begin
        model_32 = BlackScholesProcess(Float32(sigma))
        @test pricer(model_32, Underlying(Float32(S0), Float32(d)), ZeroRate(Float32(r)), mc, Smile(Float32(T), Float32.(K), isCall)) isa CuVector{Float32}
    end
end

@testset "dual numbers" begin
    smile = FMC.EuropeanOptionSmile(T, K, isCall)
    options = FMC.smile_to_options(smile)
    vega(payoffs) = ForwardDiff.derivative(s -> pricer(BlackScholesProcess(s), und, rfCurve, mc, payoffs), sigma)
    @test vega(smile) isa CuVector{Float64}
    @test Array(vega(smile)) ≈ vega(options) rtol = 1e-6
end
