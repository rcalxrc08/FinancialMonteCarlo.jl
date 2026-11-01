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
model = BlackScholesProcess(sigma, Underlying(S0, d));
K = [80.0, 90.0, 100.0, 110.0, 120.0];
isCall = [true, false, true, false, true];
mc = MonteCarloConfiguration(Nsim, Nstep, FMC.StandardMC(), FMC.CudaMode(1, 2));

@testset "$(nameof(Smile))" for Smile in (FMC.EuropeanOptionSmile, FMC.BinaryEuropeanOptionSmile)
    smile = Smile(T, K, isCall)
    options = FMC.smile_to_options(smile)
    @test which(FMC.payoff, Tuple{CuMatrix{Float64}, typeof(smile), typeof(rfCurve), typeof(mc), Float64}).module === Base.get_extension(FMC, :FinancialMonteCarloCUDAExt)

    @testset "the results stay on the GPU" begin
        @test pricer(model, rfCurve, mc, smile) isa CuVector{Float64}
        @test variance(model, rfCurve, mc, smile) isa CuVector{Float64}
        @test delta(model, rfCurve, mc, smile, 1e-4) isa CuVector{Float64}
        @test rho(model, rfCurve, mc, smile, 1e-4) isa CuVector{Float64}
        FMC.set_seed!(mc)
        @test payoff(simulate(model, rfCurve, mc, T), smile, rfCurve, mc) isa CuMatrix{Float64}
    end

    @testset "same results of the single options" begin
        @test Array(pricer(model, rfCurve, mc, smile)) ≈ pricer(model, rfCurve, mc, options)
        @test Array(variance(model, rfCurve, mc, smile)) ≈ variance(model, rfCurve, mc, options)
        @test all(collect.(confinter(model, rfCurve, mc, smile)) .≈ collect.(confinter(model, rfCurve, mc, options)))
        FMC.set_seed!(mc)
        S = simulate(model, rfCurve, mc, 2 * T)
        @test Array(payoff(S, smile, rfCurve, mc, 2 * T))[:, 2] ≈ Array(payoff(S, options[2], rfCurve, mc, 2 * T))
    end

    @testset "strikes and flags of other array types" begin
        @test pricer(model, rfCurve, mc, Smile(T, K, BitVector(isCall))) ≈ pricer(model, rfCurve, mc, smile)
        @test pricer(model, rfCurve, mc, Smile(T, CuArray(K), CuArray(isCall))) ≈ pricer(model, rfCurve, mc, smile)
    end

    @testset "Float32 stays Float32" begin
        model_32 = BlackScholesProcess(Float32(sigma), Underlying(Float32(S0), Float32(d)))
        @test pricer(model_32, ZeroRate(Float32(r)), mc, Smile(Float32(T), Float32.(K), isCall)) isa CuVector{Float32}
    end
end

@testset "dual numbers" begin
    smile = FMC.EuropeanOptionSmile(T, K, isCall)
    options = FMC.smile_to_options(smile)
    vega(payoffs) = ForwardDiff.derivative(s -> pricer(BlackScholesProcess(s, Underlying(S0, d)), rfCurve, mc, payoffs), sigma)
    @test vega(smile) isa CuVector{Float64}
    @test Array(vega(smile)) ≈ vega(options) rtol = 1e-6
end
