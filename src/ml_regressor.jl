#############################################################################
#Copyright (c) 2016 Charles Le Losq
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
function mlregressor(x::Array{Float64},y::Array{Float64},algorithm::AbstractString;X_test::Array{Float64}=[0.0],y_test::Array{Float64}=[0.0],rand_state=42,test_sz = 0.33,param_grid_kr = Dict("alpha"=> [1e1, 1e0, 0.5, 0.1, 5e-2, 1e-2, 5e-3, 1e-3],"gamma"=> logspace(-4, 4, 9)),param_grid_svm=Dict("C"=> [1e0, 2e0, 5e0, 1e1, 5e1, 1e2, 5e2, 1e3, 5e3, 1e4, 5e4, 1e5],"gamma"=> logspace(-4, 4, 9)),user_kernel="rbf")
    
	if size(X_test,1) == 1
		X_train, X_test, y_train, y_test = cross_validation.train_test_split(x, reshape(y,size(y,1),1), test_size=test_sz, random_state=rand_state) 
	elseif size(X_test,2) == size(x,2)
		X_train = x[:,:]
		y_train = y[:,:]
	else
		error("You tried to provide a testing dataset that has a different number of features (in columns) than the training set. Please correct this.")
	end
	
	# to avoid any problem with SciKit Learn quite annoying demands for the shape of arrays...
	y_train =reshape(y_train,size(y_train,1),1)
	y_test = reshape(y_test,size(y_test,1),1)
	
	# initialising the preprocessor scaler
	X_scaler = preprocessing[:StandardScaler]()
	Y_scaler = preprocessing[:StandardScaler]()
	
	X_scaler[:fit](X_train)
	Y_scaler[:fit](y_train)
	
	#scaling the data
	X_train_sc = X_scaler[:transform](X_train)
	y_train_sc = Y_scaler[:transform](y_train)

	X_test_sc = X_scaler[:transform](X_test)
	y_test_sc = Y_scaler[:transform](y_test)
	
	if algorithm=="KernelRidge"
		clf_kr = kernel_ridge[:KernelRidge](kernel=user_kernel, gamma=0.1)
		model = grid_search[:GridSearchCV](clf_kr,cv=5,param_grid=param_grid_kr)
	elseif algorithm=="SVM"
		clf_svm = svm[:SVR](kernel=user_kernel, gamma=0.1)
		model = grid_search[:GridSearchCV](clf_svm,cv=5,param_grid=param_grid_svm)
	elseif algorithm=="LinearRegression"
		model = linear_model[:LinearRegression]()
	end
	
	model[:fit](X_train_sc, squeeze(y_train_sc,2))
	predict_train_sc = model[:predict](X_train_sc)
	prediction_train = Y_scaler[:inverse_transform](predict_train_sc)
	predict_test_sc = model[:predict](X_test_sc)
	prediction_test = Y_scaler[:inverse_transform](predict_test_sc)
	
	return prediction_train, prediction_test, model, X_scaler, Y_scaler
end
