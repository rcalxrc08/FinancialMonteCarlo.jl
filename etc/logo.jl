# Draws etc/logo.png and etc/favicon.ico (copied into docs/src/assets), in the style of the FinancialToolbox.jl logo:
# three discs in the Julia colours, each drawn from what this package computes.
# Needs Luxor.jl, FinancialMonteCarlo.jl and FinancialToolbox.jl in the active environment: julia etc/logo.jl
using Luxor
using FinancialMonteCarlo
using FinancialToolbox: blsprice

const WIDTH, HEIGHT = 604, 552
const RING_RADIUS, DISC_RADIUS = 135, 117
const STROKE = 11
const BLACK = "#000000"

# Ring and disc colours of the Julia logo.
const GREEN = ("#389826", "#60ad51")
const RED = ("#cb3c33", "#d5635c")
const PURPLE = ("#9558b2", "#aa79c1")

# Market data shared by the three symbols.
const S0, K, r, d, T, σ = 100.0, 100.0, 0.0, 0.0, 1.0, 0.25
const RATES = ZeroRate(r)

function draw_disc(centre, colours)
    setcolor(colours[1])
    circle(centre, RING_RADIUS, :fill)
    setcolor(colours[2])
    circle(centre, DISC_RADIUS, :fill)
end

function draw_line(points)
    setcolor(BLACK)
    poly(points, :stroke, close = false)
end

# Maps the values in `values` linearly from `from` onto `onto`.
rescale(values, from, onto) = onto[1] .+ (values .- from[1]) .* (onto[2] - onto[1]) ./ (from[2] - from[1])

# Paths drawn after seeding the generator, as the metrics do, so the logo is the same at every run.
function seeded_simulate(model, Nsim, Nstep, seed)
    config = MonteCarloConfiguration(Nsim, Nstep, seed)
    FinancialMonteCarlo.set_seed!(config)
    return simulate(model, RATES, config, T)
end

# The simulation: Kou paths leaving the spot, a few bold ones over the crowd, and the histogram of where all the paths end.
function draw_paths()
    Nsim, Nstep = 2000, 10
    model = KouProcess(0.2, 1.0, 0.5, 12.0, 12.0, Underlying(S0, d))
    logS = log.(seeded_simulate(model, Nsim, Nstep, 18) ./ S0)
    # The drawn paths are the ones ending on evenly spaced quantiles, so the fan is spread over the whole disc.
    order = sortperm(logS[:, end])
    thin, bold = order[round.(Int, range(0.03, 0.97, length = 7) .* Nsim)], order[round.(Int, [0.1, 0.5, 0.9] .* Nsim)]
    range_y = (-0.5, 0.5)
    to_y(values) = rescale(values, range_y, (70, -70))
    xs = range(-90, 34, length = Nstep + 1)
    setline(0.3 * STROKE)
    foreach(i -> draw_line(Point.(xs, to_y(logS[i, :]))), thin)
    setline(0.8 * STROKE)
    foreach(i -> draw_line(Point.(xs, to_y(logS[i, :]))), bold)
    edges = range(range_y..., length = 13)
    counts = [count(x -> lo <= x < hi, logS[:, end]) for (lo, hi) in zip(edges[1:(end - 1)], edges[2:end])]
    lengths = 44 .* counts ./ maximum(counts)
    bar_y = to_y(edges)
    setcolor(BLACK)
    for (len, bottom, top) in zip(lengths, bar_y[1:(end - 1)], bar_y[2:end])
        len > 1 && box(Point(50, top + 1.5), Point(50 + len, bottom - 1.5), :fill)
    end
end

# The estimator: the price of a call against the number of paths, on a log scale, inside its 99% confidence interval,
# which closes on the Black and Scholes price.
function draw_convergence()
    sims = round.(Int, 10 .^ range(1.6, 3.6, length = 30))
    model = BlackScholesProcess(σ, Underlying(S0, d))
    option = EuropeanOption(T, K)
    config(Nsim) = MonteCarloConfiguration(Nsim, 1, 2)
    prices = [pricer(model, RATES, config(Nsim), option) for Nsim in sims]
    intervals = [confinter(model, RATES, config(Nsim), option) for Nsim in sims]
    exact = blsprice(S0, K, r, T, σ, d)
    half_width = maximum(abs.(collect(Iterators.flatten(intervals)) .- exact))
    to_y(price) = rescale(price, (exact - half_width, exact + half_width), (66, -66))
    xs = range(-88, 88, length = length(sims))
    setcolor(RED[1])
    poly([Point.(xs, to_y.(first.(intervals))); reverse(Point.(xs, to_y.(last.(intervals))))], :fill)
    setline(0.4 * STROKE)
    draw_line([Point(-96, 0), Point(96, 0)])
    setline(0.8 * STROKE)
    draw_line(Point.(xs, to_y.(prices)))
end

# The multivariate models: the final log returns, around their mean, of two Black and Scholes underlyings joined by a Gaussian copula.
function draw_copula()
    Nsim = 200
    marginal = BlackScholesProcess(σ, Underlying(S0, d))
    model = GaussianCopulaNVariateLogProcess(marginal, marginal, 0.7)
    S = seeded_simulate(model, Nsim, 1, 18)
    x, y = (log.(S[i][:, end]) .- log(S0 * exp((r - d - σ^2 / 2) * T)) for i in 1:2)
    pixels = 84 / (2.2 * σ)
    setcolor(BLACK)
    for (xi, yi) in zip(x, y)
        point = Point(pixels * xi, -pixels * yi)
        distance(point, O) < 88 && circle(point, 4, :fill)
    end
end

function draw_symbol(draw, centre)
    gsave()
    translate(centre)
    draw()
    grestore()
end

function draw_logo()
    background(1, 1, 1, 0)
    setline(STROKE)
    setlinecap("round")
    setlinejoin("round")
    top, left, right = Point(302, 145), Point(146, 415), Point(458, 415)
    for (centre, colours, draw) in ((top, GREEN, draw_paths), (left, RED, draw_convergence), (right, PURPLE, draw_copula))
        draw_disc(centre, colours)
        draw_symbol(draw, centre)
    end
end

# The logo as a PNG, centred on a square of the given side.
function square_png(side)
    path = tempname() * ".png"
    Drawing(side, side, path)
    setstrokescale(true)
    scale(side / WIDTH)
    translate(0, (WIDTH - HEIGHT) / 2)
    draw_logo()
    finish()
    return read(path)
end

# An ICO file is a directory of images, here PNGs, each with its size and position in the file.
function write_ico(path, sides)
    pngs = square_png.(sides)
    open(path, "w") do io
        write(io, UInt16(0), UInt16(1), UInt16(length(pngs)))
        offset = 6 + 16 * length(pngs)
        for (side, png) in zip(sides, pngs)
            write(io, UInt8(side), UInt8(side), UInt16(0), UInt16(1), UInt16(32), UInt32(length(png)), UInt32(offset))
            offset += length(png)
        end
        foreach(png -> write(io, png), pngs)
    end
end

Drawing(WIDTH, HEIGHT, joinpath(@__DIR__, "logo.png"))
draw_logo()
finish()
write_ico(joinpath(@__DIR__, "favicon.ico"), [16, 32, 48])
for asset in ("logo.png", "favicon.ico")
    cp(joinpath(@__DIR__, asset), joinpath(@__DIR__, "..", "docs", "src", "assets", asset), force = true)
end
