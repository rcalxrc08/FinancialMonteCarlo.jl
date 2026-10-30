module FinancialMonteCarloDualNumbersExt

using DualNumbers: Dual

Base.rtoldefault(::Type{T}) where {T <: Dual} = sqrt(eps(T))

end
