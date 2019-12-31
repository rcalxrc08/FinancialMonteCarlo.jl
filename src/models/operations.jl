import Base.|>;
import Base.+;

const MarketDataSet=Dict{String,FinancialMonteCarlo.AbstractMonteCarloProcess};
export MarketDataSet;

#Strategies Implementation

function |>(x::String,y::FinancialMonteCarlo.AbstractMonteCarloProcess)
	out=MarketDataSet( x => y );

	return out;
end

function |>(x::String,y::FinancialMonteCarlo.NDimensionalMonteCarloProcess)
	sep=findfirst("_",x);
	if(isnothing(sep))
		error("Nvariate process must follow the format: INDEX1_INDEX2");
	end
	len_=length(y.models);
	
	sep=split(x,"_");
	@assert len_==length(sep)
	out=MarketDataSet( x => y );
	for (key_,model_) in zip(sep,y.models)
		out[key_]=model_
	end

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
