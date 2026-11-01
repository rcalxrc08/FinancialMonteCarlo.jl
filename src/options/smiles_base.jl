function validate_abstract_smile(T::T1, K::T3, isCall::T4) where {T1, T2, T3 <: (AbstractArray{T2}), T4 <: (AbstractArray{Bool})}
    @assert length(K) == length(isCall) "K and isCall must have the same length"
    @assert minimum(K) > 0 "Strikes must be positive"
    @assert T > 0 "Maturity must be positive"
end

struct OptionSmile{Option, T1, T2, T3 <: AbstractArray{T2}, T4 <: AbstractArray{Bool}} <: AbstractSmile{Option}
    T::T1
    K::T3
    isCall::T4
    function OptionSmile{Option}(T::T1, K::T3, isCall::T4) where {Option, T1, T2, T3 <: AbstractArray{T2}, T4 <: AbstractArray{Bool}}
        validate_abstract_smile(T, K, isCall)
        return new{Option, T1, T2, T3, T4}(T, K, isCall)
    end
end

maturity(smile::AbstractSmile) = smile.T

function smile_to_options(smile::AbstractSmile{T}) where {T}
    return T.(smile.T, smile.K, smile.isCall)
end

# Undiscounted payouts of the smile for one value of the underlying: one value per strike
payout(ST::Number, smile::AbstractSmile{T}) where {T} = payout_strike.(ST, smile.K, smile.isCall, T)

# Discounted payouts along the paths: one row per path, one column per strike
function payoff(S::AbstractMatrix{num}, smile::AbstractSmile{T_}, rfCurve::abstractZeroRateCurve, mcBaseData::AbstractMonteCarloConfiguration, T1::num2 = maturity(smile)) where {T_, abstractZeroRateCurve <: AbstractZeroRateCurve, num <: Number, num2 <: Number}
    r = rfCurve.r
    T = smile.T
    ST = values_at_maturity(S, T, T1, mcBaseData)
    payoff_of_strike(K, isCall) = payout_strike.(ST, K, isCall, T_)
    payoffs = reduce(hcat, map(payoff_of_strike, smile.K, smile.isCall))

    return payoffs * exp(-integral(r, T))
end
