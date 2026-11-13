# Figures of the home page, computed with the package at build time and saved as SVG in src/assets/figures (git-ignored),
# in a light and a dark variant: assets/figures.css shows the one of the theme of the reader. The palette and the axes
# are those of the other rcalxrc08 packages.
using FinancialMonteCarlo, Distributions, Plots, Statistics
# GR draws without a display on the CI runner
ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")
gr()

const FIGURES = mkpath(joinpath(@__DIR__, "src", "assets", "figures"))

# Categorical slots in fixed order, and the ink of the text, of the axes and of the grid, per theme.
const THEMES = (light = (series = ["#2563d9", "#0a8a6f", "#b86e0b", "#0b7fb3", "#2f8a3c", "#c2456e"], text = "#1f2328", muted = "#5b6470", axis = "#b9b4ab", grid = "#e6e3de", diverging = ["#2563d9", "#f7f6f3", "#c47a12"]), dark = (series = ["#4f8cff", "#00d4aa", "#ffb84d", "#38cfff", "#7ee787", "#ff6b9d"], text = "#f5f7fa", muted = "#a0aec0", axis = "#3a4466", grid = "#252d47", diverging = ["#4f8cff", "#151b2e", "#ffb84d"]))
const MARKERS = [:circle, :rect, :utriangle, :diamond, :dtriangle, :star5]

function base_plot(theme; kwargs...)
    return plot(; size = (760, 430), background_color = :transparent, foreground_color_text = theme.muted, foreground_color_guide = theme.muted, foreground_color_axis = theme.axis, foreground_color_border = theme.axis, foreground_color_legend = nothing, legend_background_color = :transparent, legendfontcolor = theme.text, titlefontcolor = theme.text, gridcolor = theme.grid, gridalpha = 1, gridlinewidth = 1, framestyle = :axes, fontfamily = "sans-serif", titlefontsize = 12, guidefontsize = 10, tickfontsize = 9, legendfontsize = 9, titlelocation = :left, left_margin = 4Plots.mm, bottom_margin = 4Plots.mm, kwargs...)
end

function save_themes(draw, name)
    for (theme_name, theme) in pairs(THEMES)
        savefig(draw(theme), joinpath(FIGURES, "$name-$theme_name.svg"))
    end
end

const S0, r, d = 100.0, 0.02, 0.01
const zero_rate = ZeroRate(r)
underlying() = Underlying(S0, d)
# Parameters chosen for visibly different paths, not calibrated to a market.
showcase_models() = ["Heston" => HestonProcess(0.4, 0.2, 0.0, 1.5, -0.7, 0.04, underlying()), "Kou" => KouProcess(0.15, 1.0, 0.3, 20.0, 10.0, underlying()), "Variance Gamma" => VarianceGammaProcess(0.2, -0.15, 0.2, underlying())]

# A few simulated paths of three models, one panel each.
function paths_figure(theme)
    T, steps, n = 1.0, 250, 24
    panels = map(enumerate(showcase_models())) do (j, (name, model))
        S = simulate(model, zero_rate, MonteCarloConfiguration(n, steps; seed = 7), T)
        t = range(0, T; length = steps + 1)
        p = base_plot(theme; title = name, xlabel = "time", ylabel = j == 1 ? "spot" : "", legend = false, ylims = (50, 170), left_margin = 2Plots.mm)
        for i = 1:n
            plot!(p, t, S[i, :]; color = theme.series[j], linewidth = 1, alpha = 0.55)
        end
        p
    end
    return plot(panels...; layout = (1, 3), size = (960, 380), background_color = :transparent)
end

# The price of a Black-Scholes call against the number of paths, with its 95% confidence interval, against the formula.
function convergence_figure(theme)
    σ, K, T = 0.2, 100.0, 1.0
    model = BlackScholesProcess(σ, underlying())
    option = EuropeanOption(T, K)
    d1 = (log(S0 / K) + (r - d + σ^2 / 2) * T) / (σ * sqrt(T))
    exact = S0 * exp(-d * T) * cdf(Normal(), d1) - K * exp(-r * T) * cdf(Normal(), d1 - σ * sqrt(T))
    paths = round.(Int, 10 .^ range(2, 6; length = 17))
    rows = map(paths) do n
        mc = MonteCarloConfiguration(n, 1; seed = 11)
        (pricer(model, zero_rate, mc, option), confinter(model, zero_rate, mc, option, 0.95))
    end
    price, lo, hi = first.(rows), [x[2][1] for x in rows], [x[2][2] for x in rows]
    p = base_plot(theme; title = "Black–Scholes call: Monte Carlo price and 95% confidence interval, against the formula", xlabel = "paths", ylabel = "price", xscale = :log10, size = (960, 430), legend = :topright, xticks = (10 .^ (2:6), ["100", "1 000", "10 000", "100 000", "1 000 000"]))
    plot!(p, paths, lo; fillrange = hi, fillalpha = 0.22, linewidth = 0, color = theme.series[1], label = "95% confidence interval")
    plot!(p, paths, price; color = theme.series[1], linewidth = 2, marker = :circle, markersize = 4, markerstrokewidth = 0, label = "Monte Carlo")
    hline!(p, [exact]; color = theme.series[3], linewidth = 2, linestyle = :dash, label = "Black–Scholes formula")
    return p
end

# Gaussian kernel density of a sample, on a grid
function kde(x, grid; h = 1.06 * std(x) * length(x)^(-1 / 5))
    return [mean(exp.(-((g .- x) ./ h) .^ 2 ./ 2)) / (h * sqrt(2π)) for g in grid]
end

# The density of the log-return at maturity of the three models, from the simulated paths.
function distributions_figure(theme)
    T = 1.0
    grid = range(-0.9, 0.6; length = 241)
    p = base_plot(theme; title = "Density of log(S_T / S₀) at T = 1, from 20 000 simulated paths", xlabel = "log-return", ylabel = "density", size = (960, 430), legend = :topleft)
    for (j, (name, model)) in enumerate(showcase_models())
        S = simulate(model, zero_rate, MonteCarloConfiguration(20_000, 100; seed = 3), T)
        plot!(p, grid, kde(log.(S[:, end] ./ S0), grid); color = theme.series[j], linewidth = 2, label = name)
    end
    return p
end

save_themes(paths_figure, "paths")
save_themes(convergence_figure, "convergence")
save_themes(distributions_figure, "distributions")

# The hero of the home page: the density of the log-return of the Heston model as it spreads in time, from simulated
# paths (a kernel density per time step), drawn in perspective by surface.jl.
include(joinpath(@__DIR__, "surface.jl"))
let model = HestonProcess(0.4, 0.2, 0.0, 1.5, -0.7, 0.04, underlying()), T = 2.0, steps = 200
    S = simulate(model, zero_rate, MonteCarloConfiguration(20_000, steps; seed = 5), T)
    grid = collect(range(-0.8, 0.5; length = 46))
    columns = round.(Int, range(steps ÷ 8, steps; length = 34)) .+ 1
    times = (columns .- 1) .* (T / steps)
    density = reduce(vcat, [permutedims(kde(log.(S[:, c] ./ S0), grid)) for c in columns])
    hi = ceil(maximum(density))
    write_surface(FIGURES, "hero-surface", grid, times, density; zlims = (0.0, hi), xticks = [x => string(x) for x in (-0.8, -0.4, 0.0, 0.4)], yticks = [0.5 => "6M", 1.0 => "1Y", 1.5 => "18M", 2.0 => "2Y"], zticks = [z => string(round(Int, z)) for z = 1:hi], xlabel = "log-return", ylabel = "time", zlabel = "density", label = "Density of the log-return of a Heston model over time, from simulated paths")
end
