import Base.|>;
import Base.+;

#Strategies Implementation

function |>(x::String,y::FinancialMonteCarlo.AbstractMonteCarloProcess)
	out=Dict{String,FinancialMonteCarlo.AbstractMonteCarloProcess}( x => y );

	return out;
end

function |>(x::String,y::FinancialMonteCarlo.GaussianCopulaBivariateProcess)
	sep=findfirst("_",x)[1];
	if(isnothing(sep))
		error("Bivariate process must follow the format: INDEX1_INDEX2");
	end
	idx_1=x[1:(sep-1)]
	idx_2=x[(sep+1):end]
	out=Dict{String,FinancialMonteCarlo.AbstractMonteCarloProcess}( x => y );
	out[idx_1]=y.model1
	out[idx_2]=y.model2

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
