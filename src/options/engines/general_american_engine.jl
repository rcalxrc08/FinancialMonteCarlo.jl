
function payoff(S1::AbstractMatrix{num}, amPayoff::AmericanPayoff, rfCurve::abstractZeroRateCurve, mcBaseData::AbstractMonteCarloConfiguration, T1::num2 = maturity(amPayoff)) where {abstractZeroRateCurve <: AbstractZeroRateCurve, num <: Number, num2 <: Number}
    T = amPayoff.T
    NStep = mcBaseData.Nstep
    Nsim = mcBaseData.Nsim
    index1 = round(Int, T / T1 * NStep) + 1
    S = collect(S1[:, 1:index1])
    Nstep = index1 - 1
    r = rfCurve.r
    dt = T / Nstep
    # initialize vectors
    exerciseTimes = (Nstep) .* ones(Int64, Nsim)
    df_exerciseTimes = [exp(-integral(r, dt * j_)) for j_ = 0:NStep]
    #define payout
    phi(Sti::numtype_) where {numtype_ <: Number} = payout(Sti, amPayoff)
    #compute payout
    @views V = @. phi(S[:, end])
    # Backward Procedure
    s_type = eltype(S)
    b_type = typeof(zero(eltype(S)) + zero(eltype(V)))
    A_2 = Matrix{s_type}(undef, 3, 3)
    Btilde = Array{b_type}(undef, 3)
    alpha = Array{b_type}(undef, 3)
    @inbounds for j = Nstep:-1:2
        @views Tmp = @. phi(S[:, j])
        inMoneyIndexes = findall(Tmp .> 0.0)
        if !isempty(inMoneyIndexes)
            @views S_I = S[inMoneyIndexes, j]
            #- Linear Regression on Quadratic Form
            A = Matrix{eltype(A_2)}(undef, length(S_I), 3)
            @views @. A[:, 1] = 1
            @views @. A[:, 2] = S_I
            @views @. A[:, 3] = S_I^2
            b = [V[j_] * df_exerciseTimes[exerciseTimes[j_]-j+1] for j_ in inMoneyIndexes]
            A_2 .= A' * A
            LuMat = lu!(A_2)
            Btilde .= A' * b
            alpha .= LuMat \ Btilde
            #alpha=A\b;
            #Continuation Value
            CV = A * alpha
            #-- Intrinsic Value
            @views IV = Tmp[inMoneyIndexes]
            #----------
            # Find premature exercise times
            Index = findall(IV .> CV)
            @views exercisePositions = inMoneyIndexes[Index]
            # Update the outputs
            @views @. V[exercisePositions] = IV[Index]
            @views @. exerciseTimes[exercisePositions] = j - 1
        end
    end
    Out = @. V * df_exerciseTimes[exerciseTimes+1]

    return Out
end
