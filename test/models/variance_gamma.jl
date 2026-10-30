include("../testutils.jl")
using RandomNumbers, ForwardDiff

const VG_REF = (fwd = 99.09861035102614, eu = 8.368903690692187, barrier = 7.469024475794258, asian_float = 4.779663836736272)

@testset "Variance Gamma" begin
    theta1, k1 = 0.01, 0.03
    model = VarianceGammaProcess(sigma, theta1, k1, Underlying(S0, d))

    @testset "standard (Xoroshiro128Plus rng)" begin
        mc_xoro = MonteCarloConfiguration(Nsim, Nstep, FinancialMonteCarlo.SerialMode(1, RandomNumbers.Xorshifts.Xoroshiro128Plus()))
        check_prices(model, rfCurve, mc_xoro, standard_payoffs(), VG_REF; toll = 0.8)
    end

    @testset "antithetic" begin
        check_prices(model, rfCurve, mc_anti, standard_payoffs(), VG_REF; toll = 0.8)
    end

    @testset "time dependent rates and dividends" begin
        r_vec = [0.00, 0.02]
        model_curve = VarianceGammaProcess(sigma, theta1, k1, Underlying(S0, FinancialMonteCarlo.Curve([0.0, 0.02], T)))
        model_curve_3 = VarianceGammaProcess(sigma, theta1, k1, Underlying(S0, FinancialMonteCarlo.Curve([0.0, 0.01, 0.02], T)))
        implied_curve = FinancialMonteCarlo.ImpliedZeroRate(r_vec, T)
        @test occursin("VarianceGammaProcess", capture_stdout(() -> display(model_curve)))
        @test occursin("VarianceGammaProcess", capture_stdout(() -> print(model_curve)))
        @testset "forward with $(name)" for (name, m, rf, cfg) in [("scalar rate", model_curve, ZeroRate(r_vec[end]), mc), ("implied curve", model_curve, implied_curve, mc), ("zero curve", model_curve, ZeroRate(r_vec, T), mc), ("3-point dividend curve", model_curve_3, implied_curve, mc), ("antithetic", model_curve, implied_curve, mc_anti)]
            @test isfinite(pricer(m, rf, cfg, Forward(T)))
        end
        payoffs = merge(standard_payoffs(), (binary_eu = BinaryEuropeanOption(T, K),))
        check_prices(model_curve, implied_curve, mc, payoffs)
    end

    @testset "drift and subordinator law with explicit parameters" begin
        dt = T / Nstep
        @test FinancialMonteCarlo.compute_drift(model) == FinancialMonteCarlo.compute_drift(model, sigma, theta1, k1)
        @test FinancialMonteCarlo.compute_subordinator_law(model, dt) == FinancialMonteCarlo.compute_subordinator_law(model, k1, dt)
        @test FinancialMonteCarlo.compute_drift(model, Float32(sigma), Float32(theta1), Float32(k1)) isa Float32
        @test FinancialMonteCarlo.compute_subordinator_law(model, Float32(k1), Float32(dt)) isa FinancialMonteCarlo.Gamma{Float32}
        model_dual = VarianceGammaProcess(sigma, theta1, ForwardDiff.Dual(k1, 1.0), Underlying(S0, d))
        @test FinancialMonteCarlo.compute_subordinator_law(model_dual, dt) isa FinancialMonteCarlo.Gamma{<:ForwardDiff.Dual}
        @test FinancialMonteCarlo.compute_subordinator_law(model_dual, ForwardDiff.value(model_dual.κ), dt) isa FinancialMonteCarlo.Gamma{Float64}
    end

    @testset "parameter validation" begin
        @test_throws AssertionError VarianceGammaProcess(-sigma, theta1, k1, Underlying(S0, d))
        @test_throws AssertionError VarianceGammaProcess(sigma, theta1, -k1, Underlying(S0, d))
        @test_throws AssertionError VarianceGammaProcess(sigma, 10000.0, k1, Underlying(S0, d))
    end
end
