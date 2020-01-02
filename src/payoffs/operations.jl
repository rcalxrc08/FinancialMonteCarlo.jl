import Base.+;
import Base.*;
import Base.-;
import Base./;

const Position=Dict{FinancialMonteCarlo.AbstractPayoff,Number};

#Strategies Implementation
function +(x::FinancialMonteCarlo.AbstractPayoff,y::FinancialMonteCarlo.AbstractPayoff)
	out=Position( x => 1.0 );
	if haskey(out,y)
		out[x]=2.0;
	else
		out[y]=1.0;
	end
	return out;
end

function +(x::Position,y::FinancialMonteCarlo.AbstractPayoff)
	out=copy(x);
	if haskey(out,y)
		out[y]+=1.0;
	else
		out[y]=1.0;
	end
	return out;
end

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

function *(x::Position,y::Number)
	out=copy(x);
	for a in keys(out)
		out[a]*=y;
	end
	return out;
end

function *(x::FinancialMonteCarlo.AbstractPayoff,y::Number)
	return Position( x => y );
end
*(y::Number,out::FinancialMonteCarlo.AbstractPayoff)=out*y;
*(y::Number,out::Position)=out*y;
/(out::FinancialMonteCarlo.AbstractPayoff,y::Number)=out*(1.0/y);
/(out::Position,y::Number)=out*(1.0/y);
+(y::FinancialMonteCarlo.AbstractPayoff,out::Position)= out+y;
-(x::FinancialMonteCarlo.AbstractPayoff)=return -1*x;
-(x::Position)=return -1*x;
-(y::FinancialMonteCarlo.AbstractPayoff,out::FinancialMonteCarlo.AbstractPayoff)= y+(-1*out);
-(y::FinancialMonteCarlo.AbstractPayoff,out::Position)= y+(-1*out);
-(y::Position,out::FinancialMonteCarlo.AbstractPayoff)= y+(-1*out);
-(y::Position,out::Position)= y+(-1*out);

#+(x::FinancialMonteCarlo.AbstractPayoff,y::Tuple{FinancialMonteCarlo.AbstractPayoff,FinancialMonteCarlo.AbstractPayoff,typeof(+)})=(x,y,+)

#+(x::Position,FinancialMonteCarlo.AbstractPayoff)


import Base.Multimedia.display;
import Base.hash;
import Base.isequal;

hash(x::Payoff) where { Payoff <: AbstractPayoff }=sum(hash(get_parameters(x)))+hash(string(Payoff))
isequal(x::Payoff1,y::Payoff2) where { Payoff1 <: AbstractPayoff , Payoff2 <: AbstractPayoff }=hash(x)==hash(y)

function display(p::Position)
	keys_=keys(p);
	flag=0;
	for key_ in keys_
		val_=p[key_]
		iszero(val_) ? continue : nothing;
		if typeof(val_) <: Real 
			if(flag!=0)
				val_ > 0.0 ? print(+) : print(-);
			end
			if(abs(val_)!=1.0)
				flag!=0 ? print(abs(val_)) : print(val_);
				print(*);
				print(key_);
			else
				val_ > 0.0 ? nothing : print(-);
				print(key_);
			end
		else
			if(flag!=0)
				print(+)
			end
			print("(");
			print(val_);
			print(")");
			print(*);
			print(key_);
		end
		flag+=1;
	end
	println("")
end

import Base.Multimedia.print;

function print(p::Position)
	keys_=keys(p);
	flag=0;
	for key_ in keys_
		val_=p[key_]
		iszero(val_) ? continue : nothing;
		if typeof(val_) <: Real 
			if(flag!=0)
				val_ > 0.0 ? print(+) : print(-);
			end
			if(abs(val_)!=1.0)
				flag!=0 ? print(abs(val_)) : print(val_);
				print(*);
				print(key_);
			else
				val_ > 0.0 ? nothing : print(-);
				print(key_);
			end
		else
			if(flag!=0)
				print(+)
			end
			print("(");
			print(val_);
			print(")");
			print(*);
			print(key_);
		end
		flag+=1;
	end
	println("")
end