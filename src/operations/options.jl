import Base.+;
import Base.*;
import Base.-;
import Base./;

const Position=Dict{AbstractPayoff,Number};

# "sum" AbstractPayoff s
function +(x::AbstractPayoff,y::AbstractPayoff)
	out=Position( x => 1.0 );
	if haskey(out,y)
		out[x]=2.0;
	else
		out[y]=1.0;
	end
	return out;
end

# "sum" Position and AbstractPayoff
function +(x::Position,y::AbstractPayoff)
	out=copy(x);
	if haskey(out,y)
		out[y]+=1.0;
	else
		out[y]=1.0;
	end
	return out;
end

# "sum" Position s
function +(x::Position,y::Position)
	out=copy(x);
	y_keys=keys(y);
	for y_key in y_keys
		if haskey(out,y_key)
			out[y_key]+=y[y_key];
		else
			out[y_key]=y[y_key];
		end
	end
	return out;
end


# scalar multiplication
function *(x::Position,y::Number)
	out=copy(x);
	for a in keys(out)
		out[a]*=y;
	end
	return out;
end

# scalar multiplication
function *(x::AbstractPayoff,y::Number)
	return Position( x => y );
end
*(y::Number,out::AbstractPayoff)=out*y;
*(y::Number,out::Position)=out*y;
/(out::AbstractPayoff,y::Number)=out*(1.0/y);
/(out::Position,y::Number)=out*(1.0/y);
+(y::AbstractPayoff,out::Position)= out+y;
-(x::AbstractPayoff)=return -1*x;
-(x::Position)=return -1*x;
-(y::AbstractPayoff,out::AbstractPayoff)= y+(-1*out);
-(y::AbstractPayoff,out::Position)= y+(-1*out);
-(y::Position,out::AbstractPayoff)= y+(-1*out);
-(y::Position,out::Position)= y+(-1*out);

#+(x::AbstractPayoff,y::Tuple{AbstractPayoff,AbstractPayoff,typeof(+)})=(x,y,+)

#+(x::Position,AbstractPayoff)


import Base.hash;
import Base.isequal;

# Needed for Dict
hash(x::Payoff) where { Payoff <: AbstractPayoff }=sum(hash(get_parameters(x)))+hash(string(Payoff))
isequal(x::Payoff1,y::Payoff2) where { Payoff1 <: AbstractPayoff , Payoff2 <: AbstractPayoff }=hash(x)==hash(y)

# Overload for show
import Base.show;
function show(io::IO,p::Position)
	keys_=collect(keys(p));
	for idx_ in 1:length(keys_)
		key_=keys_[idx_]
		val_=p[key_]
		iszero(val_) ? continue : nothing;
		if typeof(val_) <: Real 
			if(idx_!=1)
				val_ > 0.0 ? print(io,+) : print(io,-);
			end
			if(abs(val_)!=1.0)
				idx_!=1 ? print(io,abs(val_)) : print(io,val_);
				print(io,*);
				print(io,key_);
			else
				idx_==1 && val_<0.0 ? print(io,-) : nothing;
				print(io,key_);
			end
		else
			if(idx_!=1)
				print(io,+)
			end
			print(io,"(");
			print(io,val_);
			print(io,")");
			print(io,*);
			print(io,key_);
		end
	end
	#println("")
end