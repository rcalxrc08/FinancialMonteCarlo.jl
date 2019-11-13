import Base.|>;
import Base.+;

const MarketDataSet=Dict{String,FinancialMonteCarlo.AbstractMonteCarloProcess};
export MarketDataSet;

#Strategies Implementation

function |>(x::String,y::FinancialMonteCarlo.AbstractMonteCarloProcess)
	out=MarketDataSet( x => y );

	return out;
end

function |>(x::String,y::FinancialMonteCarlo.GaussianCopulaBivariateProcess)
	sep=findfirst("_",x);
	if(isnothing(sep))
		error("Bivariate process must follow the format: INDEX1_INDEX2");
	end
	sep=sep[1];
	idx_1=x[1:(sep-1)]
	idx_2=x[(sep+1):end]
	out=MarketDataSet( x => y );
	out[idx_1]=y.model1
	out[idx_2]=y.model2

	return out;
end


function +(x::MarketDataSet,y::MarketDataSet)
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