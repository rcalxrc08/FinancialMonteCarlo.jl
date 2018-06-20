using DifferentialEquations

struct equitySpotData{T1,T2,T3<:Number}
    S0::T1
    r::T2
    d::T3
    function equitySpotData(S0::T1,r::T2,d::T3) where {T1, T2, T3 <: Number}
        if S0 <= 0.0
            error("Spot price must be positive")
        else
            return new{T1,T2,T3}(S0,r,d)
        end
    end
end 

@enum MonteCarloMode standard=1 antithetic=2

export equitySpotData,MonteCarloMode;

struct MonteCarloBaseData
	param::Dict{String,Number}
	Nsim::Integer
	Nstep::Integer
	function MonteCarloBaseData(param::Dict{String,Number},Nsim::Integer,Nstep::Integer)
        if Nsim <= 0
            error("Number of Simulations must be positive")
        elseif Nstep <= 0
            error("Number of Steps must be positive")
		else
            return new(param,Nsim,Nstep)
        end
    end
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