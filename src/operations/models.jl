import Base.+;

const MarketDataSet=Dict{String,AbstractMonteCarloProcess};
export MarketDataSet;

#Strategies Implementation

function →(x::String,y::AbstractMonteCarloProcess)
	out=MarketDataSet( x => y );

	return out;
end

function →(x::String,y::NDimensionalMonteCarloProcess)
	sep=findfirst("_",x);
	if(isnothing(sep))
		error("Nvariate process must follow the format: INDEX1_INDEX2");
	end
	len_=length(y.models);
	
	sep=split(x,"_");
	@assert len_==length(sep)
	out=MarketDataSet( x => y );

	return out;
end

function +(x::MarketDataSet,y::MarketDataSet)
	out=copy(x);
	y_keys=keys(y);
	for y_key in y_keys
		if haskey(out,y_key)|| any(out_el-> any(out_el_sub->out_el_sub==y_key,split(out_el,"_")), keys(out)) || any(y_key_el->haskey(out,y_key_el),split(y_key,"_"))
			error("check keys")
		else
			out[y_key]=y[y_key];
		end
	end
	return out;
end
