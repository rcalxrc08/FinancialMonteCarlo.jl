function gausscopulagen2!(z, Σ::AbstractMatrix{<:Real}, mcBaseData::AbstractMonteCarloConfiguration)
    rand!(host_rng(mcBaseData), MvNormal(Σ), z)
    @inbounds for i = 1:size(Σ, 2)
        d = Normal(0, sqrt(Σ[i, i]))
        @views @. z[i, :] = cdf(d, z[i, :])
    end
end