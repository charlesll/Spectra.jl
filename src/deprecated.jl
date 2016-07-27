using Base.depwarn

Base.@deprecate long(data::Array{Float64},temp::Float64,wave::Float64) tlcorrection(data::Array{Float64},temp::Float64,wave::Float64;correction="long",density=2210.0)