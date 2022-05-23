using .TaylorSeries
import Base.max,Base.min,Base.isless

value_taylor(x::Taylor1)=@inbounds @views x[0];
value_taylor(x::AbstractSeries)=@inbounds @views x[0][1];
max(x::AbstractSeries,t::Number)=ifelse(value_taylor(x)>t,x,t)
min(x::AbstractSeries,t::Number)=ifelse(value_taylor(x)<t,x,t)
max(t::Number,x::AbstractSeries)=max(x,t)
min(t::Number,x::AbstractSeries)=min(x,t)
max(x::AbstractSeries,t::AbstractSeries)=ifelse(value_taylor(x)>value_taylor(t),x,t)
min(x::AbstractSeries,t::AbstractSeries)=ifelse(value_taylor(x)<value_taylor(t),x,t)

isless(t::Number, x::AbstractSeries)=isless(t,value_taylor(x))
isless(x::AbstractSeries,t::Number)=isless(value_taylor(x),t)
isless(x::AbstractSeries,t::AbstractSeries)=isless(value_taylor(x),value_taylor(t))