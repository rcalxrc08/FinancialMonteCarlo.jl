function gausscopulagen2!(z, Σ::Matrix{Float64}, mcBaseData::AbstractMonteCarloConfiguration)
    rand!(mcBaseData.rng, MvNormal(Σ), z)
    for i = 1:size(Σ, 2)
        d = Normal(0, sqrt(Σ[i, i]))
        @views z[i, :] .= cdf.(d, z[i, :])
    end
end