include("../testutils.jl")
using DualNumbers, ForwardDiff, Statistics

@testset "Contract algebra" begin
    p = standard_payoffs()
    fwd, eu, am = p.fwd, p.eu, p.am
    put = EuropeanOption(T, K, false)
    fwd_dual = Forward(dual(T, 1.0))
    Position = FinancialMonteCarlo.Position
    model = BlackScholesProcess(sigma, Underlying(S0, d))

    @testset "positions" begin
        pos = 2.0 * (EuropeanOption(T, K) * 2.0 + fwd - fwd / 2.0 - fwd) * 2.0
        @test pos[eu] == 8.0
        @test pos[fwd] == -2.0
        @test pos[am] == 0 # a payoff that is not in the position has weight zero
        @test (-fwd)[fwd] == -1
        @test (fwd - eu)[eu] == -1
        @test (1.0 * fwd)[fwd] == (fwd * 1.0)[fwd] == ((fwd / 1.0) * 1.0)[fwd]
        @test ((2.0 * fwd) / 4.0)[fwd] == 0.5
        @test (-(1.0 * fwd))[fwd] == -1
        @test length(fwd - 1.0 * fwd_dual) == 2 # dual maturity is a different contract
        @test Position(eu => 2, put => -1) == 2eu - put
        @test eu + put == put + eu # the order of the legs does not matter for equality
        @test sum([eu, 2put, eu - eu]) == eu + 2put
        @test isempty(zero(Position)) && iszero(eu - eu)

        # `show` for positions, with real and dual weights
        for x in (pos, -fwd, fwd_dual + eu, -fwd_dual + eu, fwd_dual - eu, -fwd_dual - eu, fwd_dual * dual(1.0, 1.0) + eu * 2.0 + am, -fwd_dual * dual(1.0, 1.0) + eu * 2.0 + am, eu - eu)
            @test !isempty(sprint(show, x))
        end
        @test sprint(show, eu - eu) == "0"
        # payoff `display`/`print` overloads
        for x in (fwd, fwd_dual)
            @test occursin("Forward", capture_stdout(() -> display(x)))
            @test occursin("Forward", capture_stdout(() -> print(x)))
        end
    end

    # Regression tests of #36
    @testset "payoffs as keys: hash and equality" begin
        # measured inside a function: at the top level Julia 1.10 counts the boxing of the result
        hash_allocations(x) = @allocated hash(x)
        hash_allocations(eu)
        @test hash_allocations(eu) == 0
        # the two-argument hash, used by containers, agrees with isequal, also with an array field
        b1, b2 = BermudanOption([0.5, 1.0], K), BermudanOption([0.5, 1.0], K)
        @test isequal(b1, b2) && b1 == b2
        @test hash(b1, UInt(1)) == hash(b2, UInt(1))
        @test hash((b1,)) == hash((b2,))
        @test length(unique([b1, b2])) == 1
        # the same contract with an integer strike is the same leg
        @test length(EuropeanOption(T, 100.0) + EuropeanOption(T, 100)) == 1
        @test (EuropeanOption(T, 100.0) + EuropeanOption(T, 100))[eu] == 2
        # different contracts are different
        @test !isequal(eu, put) && hash(eu) != hash(put)
        @test !isequal(EuropeanOption(T, K), BinaryEuropeanOption(T, K))
    end

    @testset "netting" begin
        @test length(eu - eu) == 0
        @test length(0 * eu) == 0 && length(0.0 * (eu + put)) == 0
        # a netted leg no longer sets the maturity of the simulation
        long_dated = EuropeanOption(5.0, K)
        @test FinancialMonteCarlo.maturity(long_dated - long_dated + eu) == T
        @test pricer(model, rfCurve, mc_small, long_dated - long_dated + eu) == pricer(model, rfCurve, mc_small, eu)
        @test pricer(model, rfCurve, mc_small, eu - eu) == 0
        # a weight whose value is zero but that carries a derivative is kept
        @test length(eu * dual(0.0, 1.0)) == 1
        @test length(eu * ForwardDiff.Dual(0.0, 1.0)) == 1
        @test ForwardDiff.derivative(w -> pricer(model, rfCurve, mc_small, w * eu + put), 0.0) ≈ pricer(model, rfCurve, mc_small, eu)
    end

    @testset "weights, order and inference" begin
        @test FinancialMonteCarlo.Position(eu => 2, put => 1.0) isa Position{<:EuropeanOption, Float64} # promoted weights
        @test collect(keys(fwd + eu + put)) == [fwd, eu, put] # construction order
        @test (@inferred pricer(model, rfCurve, mc_small, 2eu - put)) isa Float64
        @test pricer(model, rfCurve, mc_small, 2eu - put + fwd) ≈ 2 * pricer(model, rfCurve, mc_small, eu) - pricer(model, rfCurve, mc_small, put) + pricer(model, rfCurve, mc_small, fwd)
    end

    @testset "variance and confidence interval of a position" begin
        position = 2eu - put
        FinancialMonteCarlo.set_seed!(mc_small)
        S = simulate(model, rfCurve, mc_small, T)
        path_values = 2 .* payoff(S, eu, rfCurve, mc_small) .- payoff(S, put, rfCurve, mc_small)
        @test variance(model, rfCurve, mc_small, position) ≈ var(path_values)
        lower, upper = confinter(model, rfCurve, mc_small, position)
        @test lower < pricer(model, rfCurve, mc_small, position) < upper
        @test pricer(model, rfCurve, mc_small, position) ≈ mean(path_values)
        @test FinancialMonteCarlo.delta(model, rfCurve, mc_small, position) ≈ 2 * FinancialMonteCarlo.delta(model, rfCurve, mc_small, eu) - FinancialMonteCarlo.delta(model, rfCurve, mc_small, put) atol = 1e-6
        @test rho(model, rfCurve, mc_small, position) isa Float64
    end

    @testset "portfolio on a single name" begin
        name = "ENJ"
        contracts = [p.fwd, p.eu, p.am, p.barrier, p.asian_float, p.asian_fixed, Spot()]
        portfolio = name → contracts[1]
        for c in contracts[2:end]
            portfolio += name → c
        end
        @test occursin(name, sprint(show, MIME"text/plain"(), portfolio)) # show, not an overload of display (#36)
        @test occursin(name, sprint(show, MIME"text/plain"(), name → model))
        @test collect(keys(portfolio[name])) == contracts
        @test abs(pricer(name → model, rfCurve, mc_small, portfolio) - sum(pricer(model, rfCurve, mc_small, contracts))) < 1e-8
        # positions on names that are not in the market data set are not priced
        @test pricer(name → model, rfCurve, mc_small, portfolio + ("OTHER" → eu)) == pricer(name → model, rfCurve, mc_small, portfolio)
    end
end
