
function payoff(S::AbstractMatrix{num}, payoff_::PathDependentPayoff, rfCurve::abstractZeroRateCurve, T1::num2 = maturity(payoff_)) where {abstractZeroRateCurve <: AbstractZeroRateCurve, num <: Number, num2 <: Number}
    r = rfCurve.r
    T = payoff_.T
    (Nsim, NStep) = size(S)
    NStep -= 1
    index1 = round(Int, T / T1 * NStep) + 1

    @inbounds f(S::abstractArray) where {abstractArray <: AbstractArray{num_}} where {num_ <: Number} = payout(S, payoff_)
    @inbounds payoff2 = [f(view(S, i, 1:index1)) for i = 1:Nsim]

    return payoff2 * exp(-integral(r, T))
end
