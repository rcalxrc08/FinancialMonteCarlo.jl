import Base.|>;
import Base.+;

#Strategies Implementation

function |>(x::String,y::Dict{FinancialMonteCarlo.AbstractPayoff,Number})
	out=Dict{String,Dict{FinancialMonteCarlo.AbstractPayoff,Number}}( x => y );

	return out;
end

function |>(x::String,y::FinancialMonteCarlo.AbstractPayoff)

	return x|>(1.0*y);
end


function +(x::Dict{String,Dict{FinancialMonteCarlo.AbstractPayoff,Number}},y::Dict{String,Dict{FinancialMonteCarlo.AbstractPayoff,Number}})
	out=copy(x);
	y_keys=keys(y);
	for y_key in y_keys
		if haskey(out,y_key)
			out[y_key]=out[y_key]+y[y_key];
		else
			out[y_key]=y[y_key];
		end
	end
	return out;
end
