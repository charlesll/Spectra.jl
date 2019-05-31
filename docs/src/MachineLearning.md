# Machine Learning

Spectra calls the rampy.mlregressor and rampy.mlexplorer functions to provide easy-to-use access to some machine learning algorithms from the SciKit Learn python library.

More advanced ML treatments can be done within the Julia ecosystem (e.g. using Flux, MXNet, Tensorflow etc.), but this is outsite Spectra.

```@docs
mlregressor(x::Array{Float64},y::Array{Float64};X_test::Array{Float64}=[0.0],y_test::Array{Float64}=[0.0])
mlexplorer(x::Array{Float64})
```
