import Future.randjump

abstract type AbstractZeroRateCurve end

struct SimResult{T2 <: Number}
    r::T2
    function ZeroRate(r::T2) where {T2 <: Number}
       return new{T2}(r)
    end
end
