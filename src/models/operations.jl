import Base.+;
import Base.*;
import Base.-;

#Strategies Implementation

function +(x::FinancialMonteCarlo.AbstractMonteCarloProcess,y::FinancialMonteCarlo.AbstractMonteCarloProcess)
	#out=Dict{String,FinancialMonteCarlo.AbstractMonteCarloProcess}( x => 1.0 );
	out=Dict{String,FinancialMonteCarlo.AbstractMonteCarloProcess}( x.underlying.name => x );
	if haskey(out,y.underlying.name)
		error("check keys")
	else
		out[y.underlying.name]=y;
	end
	return out;
end

function +(x::Dict{String,FinancialMonteCarlo.AbstractMonteCarloProcess},y::FinancialMonteCarlo.AbstractMonteCarloProcess)
	out=copy(x);
	if haskey(out,y.underlying.name)
		error("check keys")
	else
		out[y.underlying.name]=y;
	end
	return out;
end

function +(x::Dict{String,FinancialMonteCarlo.AbstractMonteCarloProcess},y::Dict{String,FinancialMonteCarlo.AbstractMonteCarloProcess})
	out=copy(x);
	y_keys=keys(y);
	for y_key in y_keys
		if haskey(out,y_key)
			error("check keys")
		else
			out[y_key]=y[y_key];
		end
	end
	return out;
end

+(y::FinancialMonteCarlo.AbstractMonteCarloProcess,out::Dict{String,FinancialMonteCarlo.AbstractMonteCarloProcess})= out+y;