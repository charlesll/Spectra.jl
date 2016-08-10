using Base.depwarn

Base.@deprecate long(data::Array{Float64},temp::Float64,wave::Float64) tlcorrection(data::Array{Float64},temp::Float64,wave::Float64;correction::AbstractString="long";density::Float64=2210)