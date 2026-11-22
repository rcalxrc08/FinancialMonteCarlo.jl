# FMC #85: building a model and an option inside the differentiated function (FinancialCalibrator rebuilds them per
# iteration). Run in an environment with FinancialMonteCarlo, Zygote, ForwardDiff and BenchmarkTools: julia bench_validation_85.jl <label>
using FinancialMonteCarlo, Zygote, ForwardDiff, BenchmarkTools
const label = ARGS[1]
und = Underlying(100.0, 0.01)
heston(x) = HestonProcess(x[1], x[2], x[3], x[4], x[5], x[6], und)
fh(x) = (m = heston(x); o = EuropeanOption(x[7], 100.0); m.σ * m.σ₀ + m.κ * m.θ + m.ρ * m.λ + o.T)
fk(x) = (m = KouProcess(x[1], x[2], x[3], x[4], x[5], und); m.σ + m.λ * m.p + m.λ₊ * m.λ₋)
fs(x) = (m = FinancialMonteCarlo.set_parameters!(LogNormalMixture((0.2, 0.25, 0.3), (0.3, 0.3), und), x); sum(m.η) + sum(m.λ))
xh = [0.2, 0.2, 0.0, 1.5, -0.7, 0.04, 1.0]; xk = [0.2, 3.0, 0.4, 10.0, 5.0]; xs = [0.2, 0.25, 0.3, 0.3, 0.3]
for (name, f, x) in (("Heston+option", fh, xh), ("Kou", fk, xk), ("LNM set_parameters!", fs, xs))
    tz = @belapsed Zygote.gradient($f, $x)
    tf = @belapsed ForwardDiff.gradient($f, $x)
    println(rpad(label, 7), rpad(name, 22), " Zygote ", round(tz * 1e6, digits = 2), " μs   ForwardDiff ", round(tf * 1e6, digits = 3), " μs")
end
fe(x) = (b = DoubleBarrierOption(x[1], x[2], x[3], x[4]); o = EuropeanOption(x[1], x[2]); FinancialMonteCarlo.payout(110.0, o) + b.U)
xe = [1.0, 100.0, 80.0, 120.0]
println(rpad(label, 7), rpad("payout (scalar)", 22), " Zygote ", round(@belapsed(Zygote.gradient($fe, $xe)) * 1e6, digits = 2), " μs")
fv(x) = sum(FinancialMonteCarlo.payout_untyped.(x, Ref(EuropeanOption(1.0, 100.0))))
xv = collect(range(0.5, 1.5, length = 1000))
println(rpad(label, 7), rpad("payout broadcast 1000", 22), " Zygote ", round(@belapsed(Zygote.gradient($fv, $xv)) * 1e6, digits = 2), " μs")
