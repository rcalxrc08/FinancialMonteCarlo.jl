
# Longstaff-Schwartz regression of the discounted values V on 1, S, S^2 over the paths that are in
# the money at column j, through the 3x3 normal equations (accumulated without building the design
# matrix). `discount(i)` discounts V[i] from its exercise time to the time of column j.
# Returns the coefficients, or nothing when they are not identified (fewer than 3 distinct prices
# in the money), in which case there is no early exercise at this date.
function lsm_coefficients(S::AbstractMatrix, V::AbstractVector, phi::Function, discount::Function, j::Integer)
    s_0 = zero(eltype(S))
    y_0 = zero(eltype(V)) * discount(1) * s_0
    n, s1, s2, s3, s4 = 0, s_0, s_0, s_0, s_0
    y0, y1, y2 = y_0, y_0, y_0
    @inbounds for i in axes(S, 1)
        s = S[i, j]
        phi(s) > 0 || continue
        y = V[i] * discount(i)
        s_2 = s * s
        n += 1
        s1 += s
        s2 += s_2
        s3 += s_2 * s
        s4 += s_2 * s_2
        y0 += y
        y1 += y * s
        y2 += y * s_2
    end
    n < 3 && return nothing
    A_2 = [n s1 s2; s1 s2 s3; s2 s3 s4]
    LuMat = lu!(A_2, check = false)
    issuccess(LuMat) || return nothing
    return LuMat \ [y0, y1, y2]
end

# Exercise at column j the paths whose intrinsic value is positive and above the continuation value.
function lsm_exercise!(V::AbstractVector, exerciseTimes::AbstractVector, S::AbstractMatrix, phi::Function, alpha, j::Integer)
    @inbounds for i in axes(S, 1)
        s = S[i, j]
        IV = phi(s)
        if IV > 0 && IV > alpha[1] + alpha[2] * s + alpha[3] * s^2
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
    # Backward Procedure
    for j in exercise_columns(exPayoff, Nstep, dt)
        # discount from the time of column j, (j - 1) * dt, to the exercise time
        alpha = lsm_coefficients(S, V, phi, i -> df_exerciseTimes[exerciseTimes[i]-j+2], j)
        isnothing(alpha) || lsm_exercise!(V, exerciseTimes, S, phi, alpha, j)
    end
    Out = @. V * df_exerciseTimes[exerciseTimes+1]

    return Out
end
