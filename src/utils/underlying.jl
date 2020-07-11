
abstract type AbstractUnderlying end

mutable struct UnderlyingScalar{num <: Number, num2 <: Number} <: AbstractUnderlying
	S0::num
	d::num2
	function UnderlyingScalar(S0::num_,d::num_2=0.0) where {num_ <: Number, num_2 <: Number}
		if(S0<zero(num_))
			error("Underlying starting value must be positive");
		else
			return new{num_,num_2}(S0,d)
		end
	end
end

mutable struct UnderlyingVec{num <: Number, num2 <: Number, num3 <: Number} <: AbstractUnderlying
	S0::num
	d::Curve{num2,num3}
	function UnderlyingVec(S0::num_,d::Curve{num2,num3}) where {num_ <: Number, num2 <: Number, num3 <: Number}
		if(S0<zero(num_))
			error("Underlying starting value must be positive");
		else
			return new{num_,num2,num3}(S0,d)
		end
	end
end

function Underlying(S0::num_,d::Curve{num2,num3}) where {num_ <: Number, num2 <: Number, num3 <: Number}
	return UnderlyingVec(S0,d)
end

function Underlying(S0::num_,d::num_2=0.0) where {num_ <: Number, num_2 <: Number}
	return UnderlyingScalar(S0,d)
end


export Underlying;
