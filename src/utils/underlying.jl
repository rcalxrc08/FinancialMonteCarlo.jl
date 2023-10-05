###########Underlying Types
#Abstract
abstract type AbstractUnderlying end

#Scalar dividend
struct UnderlyingScalar{num <: Number, num2 <: Number} <: AbstractUnderlying
    S0::num
    d::num2
    function UnderlyingScalar(S0::num_, d::num_2 = 0.0) where {num_ <: Number, num_2 <: Number}
        @assert S0 >= zero(num_) "Underlying starting value must be positive"
        return new{num_, num_2}(S0, d)
    end
end
#Curve dividend
struct UnderlyingVec{num <: Number, num2 <: Number, num3 <: Number} <: AbstractUnderlying
    S0::num
    d::CurveType{num2, num3}
    function UnderlyingVec(S0::num_, d::CurveType{num2, num3}) where {num_ <: Number, num2 <: Number, num3 <: Number}
        @assert S0 >= zero(num_) "Underlying starting value must be positive"
        return new{num_, num2, num3}(S0, d)
    end
end

#Common Constructors
function Underlying(S0::num_, d::CurveType{num2, num3}) where {num_ <: Number, num2 <: Number, num3 <: Number}
    return UnderlyingVec(S0, d)
end
function Underlying(S0::num_, d::num_2 = 0.0) where {num_ <: Number, num_2 <: Number}
    return UnderlyingScalar(S0, d)
end

export Underlying;
