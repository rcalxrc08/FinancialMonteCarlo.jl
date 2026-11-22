# Longstaff-Schwartz regression of the discounted values V on 1, x, x^2, with x the price relative to the spot
# (`lsm_regressor`), over the paths that are in the money at column j, through the 3x3 normal equations (accumulated
# without building the design matrix). `discount(i)` discounts V[i] from its exercise time to the time of column j.
# On the prices themselves (around 100) the columns 1, S, S^2 are almost collinear and the normal equations hold sums
# of S^4 next to the number of paths: the fitted continuation values lost two digits more than on x (#35). The fitted
# values do not depend on the change of variable, only their rounding does.
# Returns the coefficients, or nothing when they are not identified (fewer than 3 distinct prices
# in the money), in which case there is no early exercise at this date.
function lsm_coefficients(S::AbstractMatrix, V::AbstractVector, phi::Function, discount::Function, j::Integer, scale::Number)
    x_0 = zero(eltype(S)) / scale
    y_0 = zero(eltype(V)) * discount(1) * x_0
    n, x1, x2, x3, x4 = 0, x_0, x_0, x_0, x_0
    y0, y1, y2 = y_0, y_0, y_0
    @inbounds for i in axes(S, 1)
        s = S[i, j]
        phi(s) > 0 || continue
        x = lsm_regressor(s, scale)
        y = V[i] * discount(i)
        x_2 = x * x
        n += 1
        x1 += x
        x2 += x_2
        x3 = muladd(x_2, x, x3)
        x4 = muladd(x_2, x_2, x4)
        y0 += y
        y1 = muladd(y, x, y1)
        y2 = muladd(y, x_2, y2)
    end
    n < 3 && return nothing
    return solve_3x3((n * one(x_0), x1, x2, y0), (x1, x2, x3, y1), (x2, x3, x4, y2))
end

# The spot, the price of every path at t = 0, or one if it is zero (an engine starting at the origin)
function lsm_scale(S::AbstractMatrix)
    s = abs(S[1, 1])
    return iszero(s) ? one(s) : s
end
# The regressor of the price s: its return from the spot
lsm_regressor(s, scale) = (s - scale) / scale
# The continuation value fitted at the price s (Horner, with muladd)
lsm_continuation(alpha, s, scale) = evalpoly(lsm_regressor(s, scale), alpha)

# The rows with the pivot first: the row with the largest entry in column k in absolute value (the first one on a tie,
# as LAPACK) swapped with the first row
pivot_rows(k, r1, r2) = abs(r2[k]) > abs(r1[k]) ? (r2, r1) : (r1, r2)
function pivot_rows(k, r1, r2, r3)
    p1, p2 = pivot_rows(k, r1, r2)
    return abs(r3[k]) > abs(p1[k]) ? (r3, r2, r1) : (p1, p2, r3)
end
# Row r minus the multiple of the pivot row p that zeroes its column k (the columns up to k are not read afterwards)
eliminate(r, p, k) = (l = r[k] / p[k]; map((r_c, p_c) -> muladd(-l, p_c, r_c), r, p))

# The solution of the 3x3 system whose rows, augmented with the right hand side, are r1, r2, r3: Gaussian elimination
# with partial pivoting, as the LAPACK LU, on scalars, so it allocates nothing and every AD backend goes through it
# (Enzyme has no rule for the LAPACK LU). Returns nothing when the matrix is singular.
function solve_3x3(r1::Tuple{Number, Number, Number, Number}, r2::Tuple{Number, Number, Number, Number}, r3::Tuple{Number, Number, Number, Number})
    r1, r2, r3 = pivot_rows(1, r1, r2, r3)
    iszero(r1[1]) && return nothing
    r2, r3 = pivot_rows(2, eliminate(r2, r1, 1), eliminate(r3, r1, 1))
    iszero(r2[2]) && return nothing
    r3 = eliminate(r3, r2, 2)
    iszero(r3[3]) && return nothing
    x3 = r3[4] / r3[3]
    x2 = muladd(-r2[3], x3, r2[4]) / r2[2]
    x1 = muladd(-r1[3], x3, muladd(-r1[2], x2, r1[4])) / r1[1]
    return (x1, x2, x3)
end

# Exercise at column j the paths whose intrinsic value is positive and above the continuation value.
function lsm_exercise!(V::AbstractVector, exerciseTimes::AbstractVector, S::AbstractMatrix, phi::Function, alpha, j::Integer, scale::Number)
    @inbounds for i in axes(S, 1)
        s = S[i, j]
        IV = phi(s)
        if IV > 0 && IV > lsm_continuation(alpha, s, scale)
            V[i] = IV
            exerciseTimes[i] = j - 1
        end
    end
end

# The columns of the paths where the option can be exercised before maturity (column 1 is t = 0), latest first
exercise_columns(::AmericanPayoff, Nstep::Integer, ::Number) = Nstep:-1:2
exercise_columns(bmPayoff::BermudanPayoff, Nstep::Integer, dt::Number) = reverse(filter(j -> j <= Nstep, round.(Int, bmPayoff.T_ex ./ dt) .+ 1))

function payoff(S1::AbstractMatrix{num}, exPayoff::Union{AmericanPayoff, BermudanPayoff}, rfCurve::abstractZeroRateCurve, mcBaseData::AbstractMonteCarloConfiguration, T1::num2 = maturity(exPayoff)) where {abstractZeroRateCurve <: AbstractZeroRateCurve, num <: Number, num2 <: Number}
    T = maturity(exPayoff)
    index1 = maturity_index(T, T1, mcBaseData)
    S = view(S1, :, 1:index1)
    Nstep = index1 - 1
    r = rfCurve.r
    dt = T / Nstep
    # initialize vectors
    exerciseTimes = fill(Nstep, mcBaseData.Nsim)
    df_exerciseTimes = [exp(-integral(r, dt * j_)) for j_ = 0:Nstep]
    #define payout
    phi(Sti::numtype_) where {numtype_ <: Number} = payout(Sti, exPayoff)
    #compute payout
    @views V = @. phi(S[:, end])
    scale = lsm_scale(S)
    # Backward Procedure
    for j in exercise_columns(exPayoff, Nstep, dt)
        # discount from the time of column j, (j - 1) * dt, to the exercise time
        alpha = lsm_coefficients(S, V, phi, i -> df_exerciseTimes[exerciseTimes[i]-j+2], j, scale)
        isnothing(alpha) || lsm_exercise!(V, exerciseTimes, S, phi, alpha, j, scale)
    end
    Out = @. V * df_exerciseTimes[exerciseTimes+1]

    return Out
end
