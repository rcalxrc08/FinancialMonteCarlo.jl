using .DualNumbers
Base.rtoldefault(::Type{T}) where {T<:Dual} = sqrt(eps(T))