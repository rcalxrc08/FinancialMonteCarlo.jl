include("../testutils.jl")
using DualNumbers

@testset "Contract algebra" begin
    p = standard_payoffs()
    fwd, eu, am = p.fwd, p.eu, p.am
    fwd_dual = Forward(dual(T, 1.0))

    @testset "positions" begin
        pos = 2.0 * (EuropeanOption(T, K) * 2.0 + fwd - fwd / 2.0 - fwd) * 2.0
        @test pos[eu] == 8.0
        @test pos[fwd] == -2.0
        @test (-fwd)[fwd] == -1
        @test (fwd - eu)[eu] == -1
        @test (1.0 * fwd)[fwd] == (fwd * 1.0)[fwd] == ((fwd / 1.0) * 1.0)[fwd]
        @test ((2.0 * fwd) / 4.0)[fwd] == 0.5
        @test (-(1.0 * fwd))[fwd] == -1
        @test length(fwd - 1.0 * fwd_dual) == 2 # dual maturity is a different contract

        # `show` for positions, with real and dual weights
        for x in (pos, -fwd, fwd_dual + eu, -fwd_dual + eu, fwd_dual - eu, -fwd_dual - eu, fwd_dual * dual(1.0, 1.0) + eu * 2.0 + am, -fwd_dual * dual(1.0, 1.0) + eu * 2.0 + am)
            @test !isempty(sprint(show, x))
        end
        # payoff `display`/`print` overloads
        for x in (fwd, fwd_dual)
            @test occursin("Forward", capture_stdout(() -> display(x)))
            @test occursin("Forward", capture_stdout(() -> print(x)))
        end
    end

    @testset "portfolio on a single name" begin
        model = BlackScholesProcess(sigma, Underlying(S0, d))
        name = "ENJ"
        contracts = [p.fwd, p.eu, p.am, p.barrier, p.asian_float, p.asian_fixed, Spot()]
        portfolio = name → contracts[1]
        for c in contracts[2:end]
            portfolio += name → c
        end
        @test occursin(name, capture_stdout(() -> display(portfolio)))
        @test abs(pricer(name → model, rfCurve, mc, portfolio) - sum(pricer(model, rfCurve, mc, contracts))) < 1e-8
    end
end
