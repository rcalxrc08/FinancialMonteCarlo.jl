import Base.|>;
import Base.+;

#Strategies Implementation

const Portfolio=Dict{String,Position};

function |>(x::String,y::Position)
	out=Portfolio( x => y );

	return out;
end

function |>(x::String,y::FinancialMonteCarlo.SingleNamePayoff)
	sep=findfirst("_",x)[1];
	if(!isnothing(sep))
		error("NO UNDERSCORE ALLOWED IN SINGLE NAME OPTIONS");
	end
	return x|>(1.0*y);
end

function |>(x::String,y::FinancialMonteCarlo.BasketPayoff)
	sep=findfirst("_",x)[1];
	if(isnothing(sep))
		error("Bivariate payoff underlying must follow the format: INDEX1_INDEX2");
	end
	return x|>(1.0*y);
end


function +(x::Portfolio,y::Portfolio)
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

import Base.Multimedia.display;

function display(p::Portfolio)
	keys_=keys(p);
	for key_ in keys_
		print(key_);
		print(" => ")
		print(p[key_])
	end
	println("")
end
