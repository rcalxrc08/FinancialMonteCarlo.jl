function validate_abstract_smile(T::T1, K::T3, isCall::T4) where {T1, T2, T3<:(AbstractVector{T2}), T4<:(AbstractVector{Bool})}
    check_dimension(length(K) == length(isCall), "K and isCall must have the same length")
    check_domain(minimum(K) > 0, K, "Strikes must be positive")
    check_domain(T > 0, T, "Maturity must be positive")
end

struct OptionSmile{Option, T1, T2, T3 <: AbstractVector{T2}, T4 <: AbstractVector{Bool}} <: AbstractSmile{Option}
    T::T1
    K::T3
    isCall::T4
    function OptionSmile{Option}(T::T1, K::T3, isCall::T4) where {Option, T1, T2, T3 <: AbstractVector{T2}, T4 <: AbstractVector{Bool}}
        validate_abstract_smile(T, K, isCall)
        return OptionSmile{Option}(T, K, isCall, Val(:unchecked))
    end
    # Internal: no checks, for a smile already validated whose arrays only move (Adapt.adapt_structure, #73)
    OptionSmile{Option}(T::T1, K::T3, isCall::T4, ::Val{:unchecked}) where {Option, T1, T2, T3 <: AbstractVector{T2}, T4 <: AbstractVector{Bool}} = new{Option, T1, T2, T3, T4}(T, K, isCall)
end

maturity(smile::AbstractSmile) = smile.T

# The strikes and the flags move to the device of the mode, see to_mode. Moving the arrays cannot make a valid smile
# invalid, so it is not checked again: on the GPU minimum(K) would copy back to the host at every metric call (#73).
Adapt.adapt_structure(to, smile::OptionSmile{Option}) where {Option} = OptionSmile{Option}(smile.T, Adapt.adapt(to, smile.K), Adapt.adapt(to, smile.isCall), Val(:unchecked))

function smile_to_options(smile::AbstractSmile{T}) where {T}
    return T.(smile.T, smile.K, smile.isCall)
end

# Undiscounted payouts of the smile for one value of the underlying: one value per strike
payout(ST::Number, smile::AbstractSmile{T}) where {T} = payout_strike.(ST, smile.K, smile.isCall, T)

# Discounted payouts along the paths: one row per path, one column per strike, as one broadcast of the values at maturity
# against the strikes as a row (on the GPU too: one kernel for all the strikes)
payoff(S::AbstractMatrix, smile::AbstractSmile, rfCurve::AbstractZeroRateCurve, mcBaseData::AbstractMonteCarloConfiguration, T1::Number = maturity(smile)) = smile_payoffs(S, smile, rfCurve, mcBaseData, T1)
function smile_payoffs(S::AbstractMatrix, smile::AbstractSmile{T_}, rfCurve::AbstractZeroRateCurve, mcBaseData::AbstractMonteCarloConfiguration, T1::Number) where {T_}
    T = smile.T
    ST = values_at_maturity(S, T, T1, mcBaseData)
    return payout_strike.(ST, transpose(smile.K), transpose(smile.isCall), T_) .* exp(-integral(rfCurve.r, T))
end
