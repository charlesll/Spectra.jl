# Machine Learning Regression

Spectra offers a basic access to some machine learning algorithms from the SciKit Learn python library. In addition to using them for baseline fitting, you can also use them to predict an output. For instance, if you have several spectra that you pre-processed with Spectra, you can organise them with each spectra as a row of a large array. Each column will be a channel, or also called a feature, that the machine learning algorithms will look at. From this, Spectra allows you to call

 	- the Support Vector or Kernel Ridge regression algorithms from SciKit Learn;
	
	- the Linear Regression algorithm from SciKit Learn.
	
With those algorithm, you can predict a y value (for instance, the concentration of a component that seems to affect the spectral shape) that is related to changes in your spectra. This is done with the mlregressor function. A more extensive use of the SciKit Learn algorithms can be done directly with calling SciKit Learn using a PyCall instance.

```@docs
mlregressor(x::Array{Float64},y::Array{Float64},algorithm::AbstractString;X_test::Array{Float64}=[0.0],y_test::Array{Float64}=[0.0],test_sz=0.3,scaler="MinMaxScaler",rand_state=42,param_grid_kr = Dict("alpha"=> [1e1, 1e0, 0.5, 0.1, 5e-2, 1e-2, 5e-3, 1e-3],"gamma"=> logspace(-4, 4, 9)),param_grid_svm=Dict("C"=> [1e0, 2e0, 5e0, 1e1, 5e1, 1e2, 5e2, 1e3, 5e3, 1e4, 5e4, 1e5],"gamma"=> logspace(-4, 4, 9)),user_kernel="rbf")
```
