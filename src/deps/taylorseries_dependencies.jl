using .TaylorSeries
import Base.max,Base.min,Base.isless

value_taylor(x::Taylor1)=@inbounds @views x[0];
value_taylor(x::AbstractSeries)=@inbounds @views x[0][1];
!hasmethod(max,(AbstractSeries,Number)) ? (max(x::AbstractSeries,t::Number)=ifelse(value_taylor(x)>t,x,t)) : nothing
!hasmethod(max,(Number,AbstractSeries)) ? (max(t::Number,x::AbstractSeries)=max(x,t)) : nothing
!hasmethod(max,(AbstractSeries,AbstractSeries)) ? (max(x::AbstractSeries,t::AbstractSeries)=ifelse(value_taylor(x)>value_taylor(t),x,t)) : nothing
!hasmethod(min,(AbstractSeries,Number)) ? (max(x::AbstractSeries,t::Number)=ifelse(value_taylor(x)<t,x,t)) : nothing
!hasmethod(min,(Number,AbstractSeries)) ? (max(t::Number,x::AbstractSeries)=min(x,t)) : nothing
!hasmethod(min,(AbstractSeries,AbstractSeries)) ? (max(x::AbstractSeries,t::AbstractSeries)=ifelse(value_taylor(x)<value_taylor(t),x,t)) : nothing
!hasmethod(isless,(AbstractSeries,Number)) ? (isless(x::AbstractSeries,t::Number)=isless(value_taylor(x),t)) : nothing
!hasmethod(isless,(Number,AbstractSeries)) ? (isless(t::Number, x::AbstractSeries)=isless(t,value_taylor(x))) : nothing
!hasmethod(isless,(AbstractSeries,AbstractSeries)) ? (isless(x::AbstractSeries,t::AbstractSeries)=isless(value_taylor(x),value_taylor(t))) : nothing
