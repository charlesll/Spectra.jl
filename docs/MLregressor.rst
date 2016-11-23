.. _Tips:
***************************
Machine Learning Regression
***************************

Spectra offers a basic access to some machine learning algorithms from the SciKit Learn python library. In addition to using them for baseline fitting, you can also use them to predict an output. For instance, if you have several spectra that you pre-processed with Spectra, you can organise them with each spectra as a row of a large array. Each column will be a channel, or also called a feature, that the machine learning algorithms will look at. From this, Spectra allows you to call

 	- the Support Vector or Kernel Ridge regression algorithms from SciKit Learn;
	
	- the Linear Regression algorithm from SciKit Learn.
	
With those algorithm, you can predict a y value (for instance, the concentration of a component that seems to affect the spectral shape) that is related to changes in your spectra. This is done with the mlregressor function. A more extensive use of the SciKit Learn algorithms can be done directly with calling SciKit Learn using a PyCall instance.

------------
mlregressor
------------

	function mlregressor(x,y,algorithm;X_test=[0.0],y_test=[0.0],scaler="MinMaxScaler",rand_state=42,param_grid_kr = Dict("alpha"=> [1e1, 1e0, 0.1, 1e-2, 1e-3],"gamma"=> logspace(-4, 4, 5)),param_grid_svm=Dict("C"=> [1e0, 2e0, 5e0, 1e1, 5e1, 1e2, 5e2, 1e3, 5e3, 1e4, 5e4, 1e5],"gamma"=> logspace(-4, 4, 5)))
	
INPUTS

	x: Array{Float64}, the spectra organised in rows (1 row = one spectrum). The spectra should share a common X axis;

	y::Array{Float64}, the target. Only a sigle target is possible for now.

	algorithm: AbstractString, "KernelRidge" or "SVM" or "LinearRegression"

OPTIONS

	X_test: Array{Float64}, spectra organised in rows (1 row = one spectrum) that you want to use as a testing dataset. THose spectra should not be present in the x (training) dataset. The spectra should share a common X axis;
	
	y_test: Array{Float64}, the target that you want to use as a testing dataset. Those targets should not be present in the y (training) dataset;
	
	scaler: String, the type of scaling performed. Choose between MinMaxScaler or StandardScaler, see http://scikit-learn.org/stable/modules/preprocessing.html for details. Default = "MinMaxScaler";
	
	rand_state: Float64, the random seed that is used for reproductibility of the results. Default = 42;
	
	param_grid_kr: Dictionary, containg the values of the hyperparameters that should be checked by gridsearch for the Kernel Ridge regression algorithm;
	
	param_grid_svm: Dictionary, containg the values of the hyperparameters that should be checked by gridsearch for the Support Vector regression algorithm.
	
For the last two parameters, the user is refered to the documentation of SciKit Learn. See the pages:

http://scikit-learn.org/stable/modules/kernel_ridge.html

http://scikit-learn.org/stable/modules/generated/sklearn.svm.SVR.html

OUTPUTS

	prediction_train: Array{Float64}, the predicted target values for the training y dataset;
	
	prediction_test: Array{Float64}, the predicted target values for the testing y_test dataset;
	
	model: A Scikit Learn object model, see the above link for details;
	
	X_scaler: A Scikit Learn scaler object for the x values;
	
	Y_scaler: A Scikit Learn scaler object for the y values;
	
NOTES 

	For Support Vector and Kernel Ridge regressions, mlregressor performs a cross_validation search with using 5 KFold cross validators. 

	If the results are poor with Support Vector and Kernel Ridge regressions, you will have to tune the param_grid_kr or param_grid_svm dictionnary that records the hyperparameter space to investigate during the cross validation.
	
