
function payoff(S1::AbstractMatrix{num}, amPayoff::AmericanPayoff, rfCurve::abstractZeroRateCurve, T1::num2 = maturity(amPayoff)) where {abstractZeroRateCurve <: AbstractZeroRateCurve, num <: Number, num2 <: Number}
    T = amPayoff.T
    NStep = size(S1, 2) - 1
    index1 = round(Int, T / T1 * NStep) + 1
    S = collect(S1[:, 1:index1])
    S0 = first(S)
    (Nsim, Nstep) = size(S)
    Nstep -= 1
    r = rfCurve.r
    dt = T / Nstep
    # initialize vectors
    exerciseTimes = (Nstep) .* ones(Int64, Nsim)
    df_exerciseTimes = [exp.(-integral(r, dt * j_)) for j_ = 0:NStep]
    #define payout
    phi(Sti::numtype_) where {numtype_ <: Number} = payout(Sti, amPayoff)
    #compute payout
    @views V = phi.(S[:, end])
    # Backward Procedure 
    @inbounds for j = Nstep:-1:2
        @views Tmp = phi.(S[:, j])
        inMoneyIndexes = findall(Tmp .> 0.0)
        if !isempty(inMoneyIndexes)
            S_I = view(S, inMoneyIndexes, j)
            #-- Intrinsic Value
            @views IV = Tmp[inMoneyIndexes]
            #-- Continuation Value 
            #- Linear Regression on Quadratic Form
            A = [ones(length(S_I)) S_I S_I .^ 2]
            b = [V[j_] * df_exerciseTimes[exerciseTimes[j_]-j+1] for j_ in inMoneyIndexes]
            LuMat = lu!(A' * A)
            Btilde = A' * b
            alpha = LuMat \ Btilde
            #alpha=A\b;
            #Continuation Value
            CV = A * alpha
            #----------
            # Find premature exercise times
            Index = findall(IV .> CV)
            @views exercisePositions = inMoneyIndexes[Index]
            # Update the outputs
            @views V[exercisePositions] = IV[Index]
            @views exerciseTimes[exercisePositions] .= j - 1
        end
    end
    Out = V .* df_exerciseTimes[exerciseTimes.+1]

    return Out
end
