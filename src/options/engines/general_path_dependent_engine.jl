
function payoff(S::AbstractMatrix{num}, payoff_::PathDependentPayoff, rfCurve::abstractZeroRateCurve, mcBaseData::AbstractMonteCarloConfiguration, T1::num2 = maturity(payoff_)) where {abstractZeroRateCurve <: AbstractZeroRateCurve, num <: Number, num2 <: Number}
    r = rfCurve.r
    T = payoff_.T
    Nsim = mcBaseData.Nsim
    index1 = maturity_index(T, T1, mcBaseData)

    @inbounds f(S::abstractArray) where {abstractArray <: AbstractArray{num_}} where {num_ <: Number} = payout(S, payoff_)
    @inbounds payoff2 = [f(view(S, i, 1:index1)) for i = 1:Nsim]

    return payoff2 * exp(-integral(r, T))
end

# Barrier and Asian payoffs from reductions of the rows of S, with no path indexed one by one: the minimum and maximum of
# each path for the barrier options, its mean for the Asian ones. Faster than the generic engine on every device, and the
# automatic differentiation backends see a few array operations instead of one view and one call per path.
function payoff(S::AbstractMatrix{num}, payoff_::Union{BarrierPayoff, AsianPayoff}, rfCurve::AbstractZeroRateCurve, mcBaseData::AbstractMonteCarloConfiguration, T1::num2 = maturity(payoff_)) where {num <: Number, num2 <: Number}
    return payoff_by_rows(S, payoff_, rfCurve, mcBaseData, T1)
end
# The lowest and highest value of every path: the barrier payout depends on them only through comparisons with the
# barriers, so they carry no derivative (the AD extensions mark this function so, and skip its reductions)
path_range(S_::AbstractMatrix) = (vec(minimum(S_; dims = 2)), vec(maximum(S_; dims = 2)))
# The mean of every path, as the product with the vector of weights 1 / n: a matrix-vector product, which every AD
# backend differentiates as a whole (BLAS on the CPU, CUBLAS on the GPU)
row_means(S_::AbstractMatrix) = S_ * fill!(similar(S_, size(S_, 2)), inv(size(S_, 2)))
function payoff_by_rows(S::AbstractMatrix, payoff_::BarrierPayoff, rfCurve::AbstractZeroRateCurve, mcBaseData::AbstractMonteCarloConfiguration, T1::Number = maturity(payoff_))
    T = payoff_.T
    index1 = maturity_index(T, T1, mcBaseData)
    lo, hi = path_range(view(S, :, 1:index1))
    return barrier_payout.(values_at_maturity(S, T, T1, mcBaseData), lo, hi, payoff_) .* exp(-integral(rfCurve.r, T))
end
function payoff_by_rows(S::AbstractMatrix, payoff_::AsianPayoff, rfCurve::AbstractZeroRateCurve, mcBaseData::AbstractMonteCarloConfiguration, T1::Number = maturity(payoff_))
    T = payoff_.T
    index1 = maturity_index(T, T1, mcBaseData)
    A = row_means(view(S, :, 1:index1))
    return asian_payout.(values_at_maturity(S, T, T1, mcBaseData), A, payoff_) .* exp(-integral(rfCurve.r, T))
end
