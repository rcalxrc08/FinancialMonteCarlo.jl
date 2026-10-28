
function payoff(S1::AbstractMatrix{num}, bmPayoff::BermudanPayoff, rfCurve::abstractZeroRateCurve, mcBaseData::AbstractMonteCarloConfiguration, T1::num2 = maturity(bmPayoff)) where {abstractZeroRateCurve <: AbstractZeroRateCurve, num <: Number, num2 <: Number}
    T_opt = bmPayoff.T_ex
    T = maturity(bmPayoff)
    NStep = mcBaseData.Nstep
    Nsim = mcBaseData.Nsim
    index1 = round(Int, T / T1 * NStep) + 1
    S = view(S1, :, 1:index1)
    Nstep = index1 - 1
    r = rfCurve.r
    dt = T / Nstep
    # initialize
    exerciseTimes = Nstep .* ones(Int64, Nsim)
    df_exerciseTimes = [exp(-integral(r, dt * j_)) for j_ = 0:Nstep]
    phi(Sti::numtype_) where {numtype_ <: Number} = payout(Sti, bmPayoff)
    # column of each exercise date (column 1 is t = 0); exercise at maturity is already in V
    idx_opt = filter(j -> j <= Nstep, round.(Int, T_opt ./ dt) .+ 1)
    @views V = phi.(S[:, end]) #payoff
    # Backward Procedure (the regression is shared with the American engine)
    for j in reverse(idx_opt)
        alpha = lsm_coefficients(S, V, phi, i -> df_exerciseTimes[exerciseTimes[i]-j+2], j)
        isnothing(alpha) || lsm_exercise!(V, exerciseTimes, S, phi, alpha, j)
    end
    Out = @. V * df_exerciseTimes[exerciseTimes+1]

    return Out
end
