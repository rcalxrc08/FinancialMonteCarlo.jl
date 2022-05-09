using .TaylorSeries
import Base.max,Base.min,Base.isless
max(x::AbstractSeries,t::Number)=ifelse(x[0]>t,x,t)
min(x::AbstractSeries,t::Number)=ifelse(x[0]<t,x,t)
max(t::Number,x::AbstractSeries)=max(x,t)
min(t::Number,x::AbstractSeries)=min(x,t)
max(x::AbstractSeries,t::AbstractSeries)=ifelse(x[0]>t[0],x,t)
min(x::AbstractSeries,t::AbstractSeries)=ifelse(x[0]<t[0],x,t)

isless(t::Number, x::AbstractSeries)=isless(t,x[0])
isless(x::AbstractSeries,t::Number)=isless(x[0],t)
isless(x::AbstractSeries,t::AbstractSeries)=isless(x[0],t[0])