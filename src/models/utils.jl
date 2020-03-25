quantile_exp(lam::Number,rand1::Number)::Number=-log(1-rand1)/lam;


function gausscopulagen2(t::Int, Σ::Matrix{Float64},mcBaseData::AbstractMonteCarloConfiguration)
  z = rand(mcBaseData.rng,MvNormal(Σ),t)
  for i in 1:size(Σ, 2)
    d = Normal(0, sqrt(Σ[i,i]))
    @views z[i,:].= cdf.(d, z[i,:])
  end
  return z;
end

function gausscopulagen2!(z, Σ::Matrix{Float64},mcBaseData::AbstractMonteCarloConfiguration)
  rand!(mcBaseData.rng,MvNormal(Σ),z)
  for i in 1:size(Σ, 2)
    d = Normal(0, sqrt(Σ[i,i]))
    @views z[i,:].= cdf.(d, z[i,:])
  end
  return z;
end