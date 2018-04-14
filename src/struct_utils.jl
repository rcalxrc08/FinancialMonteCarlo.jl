struct equitySpotData{T1,T2,T3<:Number}
	S0::T1
	r::T2
	d::T3
end

export equitySpotData;

struct MonteCarloBaseData
	param::Dict{String,Number}
	Nsim::Integer
	Nstep::Integer
end

export MonteCarloBaseData;