struct equitySpotData
	S0::Number
	r::Number
	d::Number
end

export equitySpotData;

struct MonteCarloBaseData
	param::Dict{String,Number}
	Nsim::Integer
	Nstep::Integer
end

export MonteCarloBaseData;