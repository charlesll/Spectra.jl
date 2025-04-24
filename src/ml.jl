#############################################################################
#Copyright (c) 2016-2025 Charles Le Losq
#
#The MIT License (MIT)
#
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the #Software without restriction, including without limitation the rights to use, copy, #modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, #and to permit persons to whom the Software is furnished to do so, subject to the #following conditions:
#
#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, #INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
#############################################################################

"""
	mlregressor(x::Array{Float64},y::Array{Float64};X_test::Array{Float64}=[0.0],y_test::Array{Float64}=[0.0],test_sz=0.3,scaler="MinMaxScaler",rand_state=42)

Uses machine learning algorithms from scikit learn to perform regression between spectra and an observed variable.
This calls the rampy.mlregressor function and creates a Python object. Any algorithm parameter can be modified in the model object.

# Example
	```julia-repl
	julia> using Spectra
	julia> model = mlregressor(X,y)
	julia> model.algorithm = "SVM"
    julia> model.user_kernel = "poly"
    julia> model.fit()
	julia> y_new = model.predict(X_new)
	```

# Docstring from rampy.mlexporer

Attributes
----------

	x : {array-like, sparse matrix}, shape = (n_samples, n_features)
		Spectra; n_features = n_frequencies.
	y : array, shape = (n_samples,)
		Returns predicted values.
	X_test : {array-like, sparse matrix}, shape = (n_samples, n_features)
		spectra organised in rows (1 row = one spectrum) that you want to use as a testing dataset. THose spectra should not be present in the x (training) dataset. The spectra should share a common X axis.
	y_test : array, shape = (n_samples,)
		the target that you want to use as a testing dataset. Those targets should not be present in the y (training) dataset.
	algorithm : String,
		"KernelRidge", "SVM", "LinearRegression", "Lasso", "ElasticNet", "NeuralNet", "BaggingNeuralNet", default = "SVM"
	scaling : Bool
		True or False. If True, data will be scaled during fitting and prediction with the requested scaler (see below),
	scaler : String
		the type of scaling performed. Choose between MinMaxScaler or StandardScaler, see [http://scikit-learn.org/stable/modules/preprocessing.html](http://scikit-learn.org/stable/modules/preprocessing.html) for details. Default = "MinMaxScaler".
	test_size : float
		the fraction of the dataset to use as a testing dataset; only used if X_test and y_test are not provided.
	rand_state : Float64
		the random seed that is used for reproductibility of the results. Default = 42.
	param_kr : Dictionary
		contain the values of the hyperparameters that should be provided to KernelRidge and GridSearch for the Kernel Ridge regression algorithm.
	param_svm : Dictionary
		containg the values of the hyperparameters that should be provided to SVM and GridSearch for the Support Vector regression algorithm.
	param_neurons : Dictionary
		contains the parameters for the Neural Network (MLPregressor model in sklearn).
		Default= dict(hidden_layer_sizes=(3,),solver = 'lbfgs',activation='relu',early_stopping=True)
	param_bagging : Dictionary
		contains the parameters for the BaggingRegressor sklearn function that uses a MLPregressor base method.
		Default= dict(n_estimators=100, max_samples=1.0, max_features=1.0, bootstrap=True,
						bootstrap_features=False, oob_score=False, warm_start=False, n_jobs=1, random_state=rand_state, verbose=0)
	prediction_train : Array{Float64}
		the predicted target values for the training y dataset.
	prediction_test : Array{Float64}
		the predicted target values for the testing y_test dataset.
	model : Scikit learn model
		A Scikit Learn object model, see scikit learn library documentation.
	X_scaler :
		A Scikit Learn scaler object for the x values.
	Y_scaler :
		A Scikit Learn scaler object for the y values.

# Remarks

For details on hyperparameters of each algorithms, please directly consult the documentation of SciKit Learn at:

[http://scikit-learn.org/stable/](http://scikit-learn.org/stable/)

For Support Vector and Kernel Ridge regressions, mlregressor performs a cross_validation search with using 5 KFold cross validators.

If the results are poor with Support Vector and Kernel Ridge regressions, you will have to tune the param_grid_kr or param_grid_svm dictionnary that records the hyperparameter space to investigate during the cross validation.

Results for machine learning algorithms can vary from run to run. A way to solve that is to fix the random_state.
For neural nets, results from multiple neural nets (bagging technique) may also generalise better, such that
it may be better to use the BaggingNeuralNet function.

"""
function mlregressor(
    x::Array{Float64},
    y::Array{Float64};
    X_test::Array{Float64}=[0.0],
    y_test::Array{Float64}=[0.0],
)
    return rampy.mlregressor(x, y; X_test=X_test, y_test=y_test)
end

"""

	mlexplorer(x::Array{Float64})

Use machine learning algorithms from scikit learn to explore spectroscopic datasets. Performs automatic scaling and train/test split before NMF or PCA fit.

# Example

	```julia-repl
	julia> explo = mlexplorer(X) # X is an array of signals built by mixing two partial components
	julia> explo.algorithm = "NMF" # using Non-Negative Matrix factorization
	julia> explo.nb_compo = 2 # number of components to use
	julia> explo.test_size = 0.3 # size of test set
	julia> explo.scaler = "MinMax" # scaler
	julia> explo.fit() # fitting!
	julia> W = explo.model.transform(explo.X_train_sc) # getting the mixture array
	julia> H = explo.X_scaler.inverse_transform(explo.model.components_) # components in the original space
	julia> plot(X,H.T) # plot the two components
	```

# Docstring from rampy.mlexporer

Attributes
----------

	x : {array-like, sparse matrix}, shape = (n_samples, n_features)
		Spectra; n_features = n_frequencies.
	X_test : {array-like, sparse matrix}, shape = (n_samples, n_features)
		spectra organised in rows (1 row = one spectrum) that you want to use as a testing dataset. THose spectra should not be present in the x (training) dataset. The spectra should share a common X axis.
	algorithm : String,
		"PCA", "NMF", default = "PCA"
	scaling : Bool
		True or False. If True, data will be scaled prior to fitting (see below),
	scaler : String
		the type of scaling performed. Choose between MinMaxScaler or StandardScaler, see [http://scikit-learn.org/stable/modules/preprocessing.html](http://scikit-learn.org/stable/modules/preprocessing.html) for details. Default = "MinMaxScaler".
	test_size : float
		the fraction of the dataset to use as a testing dataset; only used if X_test and y_test are not provided.
	rand_state : Float64
		the random seed that is used for reproductibility of the results. Default = 42.
	model : Scikit learn model
		A Scikit Learn object model, see scikit learn library documentation.

# Remarks

For details on hyperparameters of each algorithms, please directly consult the documentation of SciKit Learn at:

[http://scikit-learn.org/stable/](http://scikit-learn.org/stable/)

Results for machine learning algorithms can vary from run to run. A way to solve that is to fix the random_state.

"""
function mlexplorer(x::Array{Float64})
    return rampy.mlexplorer(x)
end
