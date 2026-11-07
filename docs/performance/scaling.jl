# Price of a European call under Black-Scholes against the number of paths: CPU, and GPU in double and single
# precision. Writes docs/src/assets/performance_scaling.svg, shown in docs/src/manual/performance.md.
# Needs a CUDA GPU and an environment with FinancialMonteCarlo, CUDA, BenchmarkTools and Plots:
#     julia --project=<that environment> -t 1 docs/performance/scaling.jl
using FinancialMonteCarlo, CUDA, BenchmarkTools, Plots, Printf
const FMC = FinancialMonteCarlo
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 2

Nstep = 50
paths = [1_000, 10_000, 100_000, 1_000_000]
function price(U, mode)
    model = BlackScholesProcess(U(2) / 10, Underlying(U(100), U(1) / 100))
    curve = ZeroRate(U(2) / 100)
    option = EuropeanOption(one(U), U(100))
    times = map(paths) do Nsim
        mc = MonteCarloConfiguration(Nsim, Nstep; mode = mode, seed = 1)
        minimum(@benchmark CUDA.@sync(pricer($model, $curve, $mc, $option)) evals = 1).time / 1e6
    end
    GC.gc()
    CUDA.reclaim()
    return times
end
series = [("CPU, Float64", "#2a78d6", price(Float64, FMC.SerialMode())), ("GPU, Float64", "#eb6834", price(Float64, FMC.CudaMode())), ("GPU, Float32", "#1baf7a", price(Float32, FMC.CudaMode()))]
for (name, _, times) in series
    println(rpad(name, 14), join([@sprintf("%10.2f ms", t) for t in times]))
end

plt = plot(; xscale = :log10, yscale = :log10, size = (720, 420), framestyle = :axes, grid = :y, gridalpha = 0.25, legend = :topleft, legendfontsize = 9, guidefontsize = 10, titlefontsize = 11, xlabel = "number of paths ($Nstep time steps)", ylabel = "time (ms)", title = "Black-Scholes European call", xlims = (paths[1] / 1.5, paths[end] * 4), xticks = (paths, ["1k", "10k", "100k", "1M"]), yticks = ([1, 10, 100, 1000], ["1", "10", "100", "1000"]))
for (name, color, times) in series
    plot!(plt, paths, times; label = name, color = color, linewidth = 2, marker = :circle, markersize = 5, markerstrokewidth = 0)
    annotate!(plt, paths[end] * 1.15, times[end], text(name, 9, :left, "#52514e"))
end
savefig(plt, joinpath(@__DIR__, "..", "src", "assets", "performance_scaling.svg"))
