# Shared fixtures and helpers for the test suite.
# Every test file includes this, so each one can also be run on its own, e.g.
#     julia --project=test test/models/black_scholes.jl
using Test, FinancialMonteCarlo

# Market data
const S0 = 100.0
const K = 100.0
const r = 0.02
const T = 1.0
const d = 0.01
const D = 90.0 # down barrier
const sigma = 0.2

# Simulation settings
const Nsim = 10_000
const Nstep = 30
const mc = MonteCarloConfiguration(Nsim, Nstep)
const mc_anti = MonteCarloConfiguration(Nsim, Nstep, FinancialMonteCarlo.AntitheticMC())
const rfCurve = ZeroRate(r)

# Black-Scholes reference prices for the market data above (also the limit of the mixture models)
const BS_CALL_REF = (fwd = 99.1078451563562, eu = 8.43005524824866, am = 8.450489415187354, barrier = 7.5008664470880735, asian_float = 4.774451704549382)
const BS_PUT_REF = (eu = 7.342077422567968, am = 7.467065889603365, barrier = 0.29071929104261723, asian_float = 4.230547012372306, asian_fixed = 4.236220218194027)

"Call payoffs on (T, K, D) priced by most model tests."
standard_payoffs() = (fwd = Forward(T), eu = EuropeanOption(T, K), am = AmericanOption(T, K), barrier = BarrierOptionDownOut(T, K, D), asian_float = AsianFloatingStrikeOption(T), asian_fixed = AsianFixedStrikeOption(T, K))

"Put counterparts of `standard_payoffs`."
standard_put_payoffs() = (eu = EuropeanOption(T, K, false), am = AmericanOption(T, K, false), barrier = BarrierOptionDownOut(T, K, D, false), asian_float = AsianFloatingStrikeOption(T, false), asian_fixed = AsianFixedStrikeOption(T, K, false))

"""
    check_prices(model, rfCurve, mc, payoffs, refs = (;); toll = 0.0)

Price every payoff in the named tuple `payoffs`, one testset each. Payoffs with an entry in
`refs` must be within `toll` of it; the others must just produce a finite price.
"""
function check_prices(model, rfCurve, mc, payoffs, refs = (;); toll = 0.0)
    for (name, payoff_) in pairs(payoffs)
        @testset "$name" begin
            price = pricer(model, rfCurve, mc, payoff_)
            if haskey(refs, name)
                @test abs(price - refs[name]) < toll
            else
                @test isfinite(price)
            end
        end
    end
end

"Return what `f()` writes to stdout (the package's `display`/`print` overloads write there directly)."
function capture_stdout(f)
    mktemp() do path, io
        redirect_stdout(f, io)
        flush(io)
        return read(path, String)
    end
end
