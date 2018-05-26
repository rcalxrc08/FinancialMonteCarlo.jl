using DifferentialEquations

struct equitySpotData{T1,T2,T3<:Number}
	S0::T1
	r::T2
	d::T3
end

@enum MonteCarloMode standard=1 antithetic=2

export equitySpotData,MonteCarloMode;

struct MonteCarloBaseData
	param::Dict{String,Number}
	Nsim::Integer
	Nstep::Integer
end

export MonteCarloBaseData;


function concrete_types(type1::Type)::Array{Type}
	subtypes_1=subtypes(type1);
	if isempty(subtypes_1)
		return [type1];
	else
		subsubtypes=Type[];
		for type1 in subtypes_1
			subtype1=concrete_types(type1)
			append!(subsubtypes,subtype1);
		end
		return subsubtypes
	end
end

export concrete_types;