
# The column of S, simulated up to T1, holding the values at maturity T
function values_at_maturity(S::AbstractMatrix, T, T1, mcBaseData::AbstractMonteCarloConfiguration)
    index1 = round(UInt, T / T1 * mcBaseData.Nstep) + 1
    return @views S[:, index1]
end

function payoff(S::AbstractMatrix{num}, euPayoff::EuropeanPayoff, rfCurve::abstractZeroRateCurve, mcBaseData::AbstractMonteCarloConfiguration, T1::num2 = maturity(euPayoff)) where {abstractZeroRateCurve <: AbstractZeroRateCurve, num <: Number, num2 <: Number}
    r = rfCurve.r
    T = euPayoff.T
    ST = values_at_maturity(S, T, T1, mcBaseData)
    phi(Sti::numtype_) where {numtype_ <: Number} = payout(Sti, euPayoff)
    payoff2 = @. phi(ST)

    return payoff2 * exp(-integral(r, T))
end
