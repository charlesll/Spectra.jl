.. _Tips:
***************************
Machine Learning Regression
***************************

Spectra offers a basic access to some machine learning algorithms from the SciKit Learn python library. In addition to using them for baseline fitting, you can also use them to predict an output. For instance, saying you have several spectra that you pre-processed with Spectra, you can organise them with each spectra as a row of a large array. Each column will be a channel, or also called a feature, that the machine learning algorithms will look at. From this, Spectra allows you to call either the Support Vector Machines or Kernel Ridge regression algorithms from SciKit Learn, in order to predict a y value (for instance, the concentration of a component that seems to affect the spectral shape). This is done with the MLregressor function.

------------
MLregressor
------------

	function mlregressor(x::Array{Float64},y::Array{Float64},algorithm::AbstractString;X_test::Array{Float64}=[0.0],y_test::Array{Float64}=[0.0],rand_state=42,test_sz = 0.33,param_grid_kr = Dict("alpha"=> [1e1, 1e0, 0.1, 1e-2, 1e-3],"gamma"=> logspace(-4, 4, 5)),param_grid_svm=Dict("C"=> [1e0, 1e1, 1e2, 1e3, 1e4],"gamma"=> logspace(-4, 4, 5)))
	
INPUTS

	x::Array{Float64}








OPTIONS

OUTPUTS


	
