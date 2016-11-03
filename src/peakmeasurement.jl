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
# functions.jl contains several mathematic functions
#
#
#############################################################################

# this function smooths a peak and measures its height and its half-width
function peakhw(x::Array{Float64},y::Array{Float64},p=0.0001;y_smo_out=false)
    ### PRELIMINARY CHECK: INCREASING SIGNAL
    if x[end,1] < x[1,1]
        x = flipdim(x,1)
    y = flipdim(y,1)
    end

    # To help the splines and the machine learning methods, we scale all the dataset between 0 and 1
  	x = reshape(x,size(x,1),1)
  	y = reshape(y,size(y,1),1)

    # First a small smoothing of the spectrum
    # Initialising the preprocessor scaler
  	X_scaler = preprocessing[:StandardScaler]()
  	Y_scaler = preprocessing[:StandardScaler]()

  	X_scaler[:fit](x)
  	Y_scaler[:fit](y)

  	# Scaling the data
  	x_sc = X_scaler[:transform](x)
  	y_sc = Y_scaler[:transform](y)

    clf = svm[:SVR](kernel="rbf", gamma=0.1)
		svr = grid_search[:GridSearchCV](clf,cv=5,param_grid=Dict("C"=> [1e-1, 1e0, 1e1, 1e2, 1e3],"gamma"=> logspace(-4, 4, 9))) # GridSearchCV for best parameters
		svr[:fit](x_sc, squeeze(y_sc,2)) #SciKit learn is expecting a y vector, not an array...
		y_smo_sc = svr[:predict](x_sc)

    y_smo = Y_scaler[:inverse_transform](y_smo_sc)

    x_maximum = x[y_smo .== maximum(y_smo)]
    x_1 = x[x .<x_maximum]
    x_2 = x[x .>=x_maximum]
    y_first_portion = y_smo[x .<x_maximum]
    y_second_portion = y_smo[x .>=x_maximum]
    half_int = maximum(y_smo)/2
    idx_1 = findmin(abs(y_first_portion-half_int))
    idx_2 = findmin(abs(y_second_portion-half_int))
    hwhm = x_2[idx_2[2]]-x_1[idx_1[2]]

    if y_smo_out == true
      return x_maximum, hwhm, y_smo
    elseif y_smo_out ==false
      return x_maximum, hwhm
    else
      error("Set y_smo_out to true or false.")
    end

end
