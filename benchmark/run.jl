# Run the benchmark suite and print a table of the results:
#     julia --project=benchmark -t auto benchmark/run.jl [group ...]
# The arguments select groups of benchmarks/benchmarks.jl (models, methods, payoffs, markets, greeks), all by default.
# BENCHMARK_SECONDS sets the time budget of each benchmark (2 s by default). On GitHub Actions the table is also written to the summary of the job.
using BenchmarkTools, Printf, Statistics
include(joinpath(@__DIR__, "benchmarks.jl"))

function run_suite(suite::BenchmarkGroup, groups::AbstractVector{<:AbstractString}; seconds::Real = 2)
    unknown = setdiff(groups, keys(suite))
    isempty(unknown) || error("unknown benchmark groups: $(join(unknown, ", ")); available: $(join(keys(suite), ", "))")
    selected = isempty(groups) ? suite : BenchmarkGroup([group => suite[group] for group in groups]...)
    for (_, benchmark) in BenchmarkTools.leaves(selected)
        benchmark.params.seconds = seconds
    end
    @info "Tuning and running $(length(BenchmarkTools.leaves(selected))) benchmarks on $(Threads.nthreads()) threads"
    tune!(selected)
    return run(selected; verbose = true)
end

format_time(ns) = ns < 1e3 ? @sprintf("%.0f ns", ns) : ns < 1e6 ? @sprintf("%.1f μs", ns / 1e3) : ns < 1e9 ? @sprintf("%.2f ms", ns / 1e6) : @sprintf("%.2f s", ns / 1e9)
format_memory(bytes) = bytes < 1024 ? "$bytes B" : bytes < 1024^2 ? @sprintf("%.1f KiB", bytes / 1024) : @sprintf("%.1f MiB", bytes / 1024^2)

# One row per benchmark, in the order of the suite: group, name, median and minimum time, memory, allocations
function results_table(results::BenchmarkGroup)
    header = ["Group", "Benchmark", "Median", "Minimum", "Memory", "Allocations"]
    rows = Vector{String}[]
    for group in sort!(collect(keys(results))), name in sort!(collect(keys(results[group])))
        trial = results[group][name]
        med, low = median(trial), minimum(trial)
        push!(rows, [group, name, format_time(time(med)), format_time(time(low)), format_memory(memory(med)), string(allocs(med))])
    end
    return header, rows
end

function markdown_table(header, rows)
    widths = [maximum(length, [header[c]; getindex.(rows, c)]) for c in eachindex(header)]
    left = [true, true, false, false, false, false]
    cell(text, c) = left[c] ? rpad(text, widths[c]) : lpad(text, widths[c])
    line(cells) = "| " * join((cell(text, c) for (c, text) in enumerate(cells)), " | ") * " |"
    rule = "|" * join((left[c] ? ":" * "-"^(widths[c] + 1) : "-"^(widths[c] + 1) * ":" for c in eachindex(header)), "|") * "|"
    return join([line(header); rule; line.(rows)], "\n")
end

results = run_suite(SUITE, ARGS; seconds = parse(Float64, get(ENV, "BENCHMARK_SECONDS", "2")))
table = markdown_table(results_table(results)...)
context = "Julia $VERSION, $(Sys.cpu_info()[1].model), $(Threads.nthreads()) threads, $Nsim paths × $Nstep steps"
println("\n", context, "\n\n", table)
summary = get(ENV, "GITHUB_STEP_SUMMARY", "")
isempty(summary) || open(io -> println(io, "## Benchmarks\n\n", context, "\n\n", table), summary, "a")
