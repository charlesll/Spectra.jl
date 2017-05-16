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
# baseline.jl contains several "baseline" functions. It is directly dependent on functions.jl.
#
#############################################################################

using JuMP
using Ipopt

"""
	baseline(x::Array{Float64},y::Array{Float64},roi::Array{Float64},basetype::AbstractString;p=1.0,SplOrder=3,roi_out="no")
	
Baseline subtraction can be made with using the baseline function:

INPUTS:

	x: Array{Float64}, containing the x values;
	
	y: Array{Float64}, containing the y values;
	
	roi: Array{Float64}, containing the "region of interest", i.e. the places where you want to fit the baseline. For instance, if the baseline should fit the regions comprised between 750 and 800 cm^{-1}, and 1250 and 1300 cm^{-1}: roi = [750. 800.; 1250. 1300.];

	basetype: AbstractString, the type of baseline that you want to use. For now, polynomial and cubic spline baselines are available. Indicate the type you want as:

		Polynomial baseline: enter "poly" for basetype, then the polynomial degree as p.

		Dierckx cubic spline baseline: enter "Dspline" for basetype, then the smoothing degree as p.

		Generalised Cross-Validated baseline: enter "gsvspline" for basetype, then the smoothing degree as p. 
		
		Kernel Ridge Regression: enter "KRregression" for basetype, no need to provide p.
		
		Support Vector Machines regression: enter "SVMregression" for basetype, no need to provide p.

OPTIONS:

	p:: Float64, if using gcvspline or Dspline, this number indicates the spline smoothing coefficient. If using "poly", it is the degree of the polynomial function to be fitted. Please enter a float number (1.0, 2.0 or 3.0 for splines of order 1, 2 or 3), and it is automatically converted to an Integer for the polyfit function. Default = 1.0.

	SplOrder: Integer, the spline coefficient to be used with the Dspline or gcvspline options. Default = 3.
	
	roi_out: String, "no" or "yes". This will result in an additional output matrix containing the y signal in the roi regions of interest, which can then be used to plot and to evaluate the roi provided to the baseline function.
	
OUTPUTS:

(are combined in a tuple in one array if only one output variable is provided)

	y_corr: Array{Float64}, the spectrum corrected from its baseline;
	
	bass: Array{Float64}, the baseline.
	
OPTIONAL OUTPUT:

	y_roi_out: Array{Float64}, an 2 column array containing the initial x-y pairs of the signal in the roi regions of interest.
	
NOTES:

Errors on measurenements are automatically provided as sqrt(y) in gcvspline. For further options, please use the gcvspl and splderivative functions that directly call the GCVSPL and SPLDER function of the gcvspl.f program (Holtring, 1986). Further informations for the use of splines are given in the Splines section, see :ref:`Splines`.

The Kernel Ridge and Support Vector Machines regression algorithms call the Scikit Learn library, available in Python. This library thus SHOULD be installed. They are machine learning algorithms that will try to automatically fit the baseline in the provided regions of interest. They are slower that splines, but have the advantage of avoiding the (sometimes painful) tuning of the spline coefficients.

The Kernel Ridge and Support Vector Machines regression algorithms used a Cross-Validated approach to increase the generalisation and avoid overfitting. The GridSearchCV function of SciKit Learn is called, with 5 fold cross-validation and the following gridsearch parameters:

	- For KRregression: param_grid=Dict("alpha"=> [1e0, 0.1, 1e-2, 1e-3],"gamma"=> logspace(-4, 4, 9));
	- For SVMregression: param_grid=Dict("C"=> [1e0, 1e1, 1e2, 1e3],"gamma"=> logspace(-4, 4, 9)).
	
Please see the SciKit Learn documentation at http://scikit-learn.org/stable/index.html for further details on the implementation of those technics, together with the source code of Spectra.jl.

EXAMPLES:

For instance, for subtracting a constant baseline between 1250 and 1300 cm^{-1}:

    roi = [1250. 1300.]
	
    basetype = "poly"
	
    y_corr, bas = baseline(x,y,roi,"poly",p=0.0)
	

For a linear baseline,
	
    bas = baseline(x,y,roi,"poly",p=1.0)

For a second order polynomial baseline,
	
    bas = baseline(x,y,roi,"poly",p=2.0)

with the last coefficient will be the one in front of x^2. This can continue as you want by adding more 1.0 values to p.

For a cubic spline baseline fitting the basis of a peak centered at 1100 cm^{-1} and with basis at 900 and 1250 cm^{-1}:

    roi = [890. 910.; 1250. 1300.]
	
    basetype = "Dspline"
	
    bas = baseline(x,y,roi,basetype,p=0.01)

p there is the smoothing parameter used. The cubic spline uses the Dierckx package initially written in Fortran and used in Julia: https://github.com/kbarbary/Dierckx.jl

"""

function baseline(x::Array{Float64},y::Array{Float64},roi::Array{Float64},basetype::AbstractString;p=1.0,SplOrder=3,roi_out="no")
    ### PRELIMINARY CHECK: INCREASING SIGNAL
	if x[end,1] < x[1,1]
	    x = flipdim(x,1)
		y = flipdim(y,1)
	end
	
	#### PRELIMINARY STEP: FIRST WE GRAB THE GOOD SIGNAL IN THE ROI
    interest_index::Array{Int64} = find(roi[1,1] .<= x[:,1] .<= roi[1,2])
    if size(roi)[1] > 1
        for i = 2:size(roi)[1]
            interest_index = vcat(interest_index,  find(roi[i,1] .<= x[:,1] .<= roi[i,2]))
        end
    end
	
    interest_x = x[interest_index,1]
    interest_y = y[interest_index,1]
	
	# To help the splines and the machine learning methods, we scale all the dataset between 0 and 1
	interest_x = reshape(interest_x,size(interest_x,1),1)
	interest_y = reshape(interest_y,size(interest_y,1),1)
	
	# Initialising the preprocessor scaler
	X_scaler = preprocessing[:StandardScaler]()
	Y_scaler = preprocessing[:StandardScaler]()

	
	X_scaler[:fit](interest_x)
	Y_scaler[:fit](interest_y)

	# Scaling the data
	x_bas_sc = X_scaler[:transform](interest_x)
	y_bas_sc = Y_scaler[:transform](interest_y)
	x_sc = X_scaler[:transform](reshape(x,size(x,1),1))
	
	# we assume errors as sqrt(y), we only use it for gcvspline
	#ese_interest_y = sqrt(abs(y_bas_sc)) 
	#ese_interest_y[interest_y.==0.0] = mean(ese_interest_y) # We protect it from any 0.0 values in y
	ese_interest_y = ones(size(y_bas_sc))
	
	#### THEN WE GO TO THE RIGHT METHOD FOR BASELINE CALCULATION

	######## POLYNOMIAL BASELINE
	if basetype == "poly"
        # The model for fitting baseline to roi signal
        #mod = Model(solver=IpoptSolver(print_level=0))
        #n::Int = size(p)[1] # number of coefficients
        #m::Int = size(interest_x)[1] # number of data

        #@variable(mod,p_val[i=1:n])
        #setvalue(p_val[i=1:n], p[i])
        #@NLobjective(mod,Min,sum{( sum{p_val[i]*interest_x[j]^(i-1), i=1:n} - interest_y[j])^2, j=1:m})
        #status = solve(mod)
        #println("Solver status: ", status)
		#best_p::Vector{Float64} = getvalue(p_val)
		best_p::Array{Float64} = polyfit(x_bas_sc[:,1],y_bas_sc[:,1],round(Int,p))
        y_calc_sc::Array{Float64} = poly(best_p,x_sc[:,1])

	######## DIERCKX SPLINE BASELINE
	elseif basetype == "Dspline"
		spl = Spline1D(x_bas_sc[:,1],y_bas_sc[:,1],s=p[1],bc="extrapolate",k=SplOrder)
		y_calc_sc = evaluate(spl,x_sc[:,1])

	######## GCV SPLINE BASELINE
	elseif basetype == "gcvspline"
		c, WK, IER = gcvspl_julia(x_bas_sc[:,1],y_bas_sc[:,1],ese_interest_y[:,1],p[1];SplineOrder = Int32(SplOrder-1)) # with cubic spline as the default
		y_calc_sc = splderivative_julia(x_sc[:,1],x_bas_sc[:,1],c,SplineOrder= Int32(SplOrder-1))

	######## KERNEL RIDGE REGRESSION WITH SCIKIT LEARN
	elseif basetype == "KRregression"
		
		clf = kernel_ridge[:KernelRidge](kernel="rbf", gamma=0.1)
		kr = model_selection[:GridSearchCV](clf,cv=5,param_grid=Dict("alpha"=> [1e1, 1e0, 0.5, 0.1, 5e-2, 1e-2, 5e-3, 1e-3,1e-4],"gamma"=> logspace(-4, 4, 9)))# GridSearchCV for best parameters
		kr[:fit](x_bas_sc, squeeze(y_bas_sc,2)) #SciKit learn is expecting a y vector, not an array...
		y_calc_sc = kr[:predict](x_sc)
		
	######## SUPPORT VECTOR MACHINES REGRESSION WITH SCIKIT LEARN
	elseif basetype == "SVMregression"
	
		clf = svm[:SVR](kernel="rbf", gamma=0.1)
		svr = model_selection[:GridSearchCV](clf,cv=5,param_grid=Dict("C"=> [1e-1, 1e0, 1e1, 1e2, 1e3],"gamma"=> logspace(-4, 4, 9))) # GridSearchCV for best parameters
		svr[:fit](x_bas_sc, squeeze(y_bas_sc,2)) #SciKit learn is expecting a y vector, not an array...
		y_calc_sc = svr[:predict](x_sc)
		
	######## GRADIENT BOOSTING ENSEMBLE REGRESSION WITH SCIKIT LEARN
	elseif basetype == "GPregression"

		# constructing a GridSearchCV instance for grabing the best parameters
		gpr = gaussian_process[:GaussianProcess](corr="squared_exponential", theta0=1e-1,thetaL=1e-3, thetaU=1,nugget=(ese_interest_y[:] ./ interest_y[:]).^2,random_start=100)
		#svr = grid_search[:GridSearchCV](clf,cv=5,param_grid=Dict("C"=> [1e-1, 1e0, 1e1, 1e2, 1e3],"gamma"=> logspace(-2, 2, 5)))
		gpr[:fit](interest_x, squeeze(interest_y,2)) #SciKit learn is expecting a y vector, not an array...
		y_calc = gpr[:predict](x_sc)
		
	######## RAISING ERROR IF NOT THE GOOD CHOICE
	else
        error("Not implemented, choose between poly, Dspline, gcvspline, KRregression and SVMregression.") # for now GPregression is hidden
    end
	
	y_calc = Y_scaler[:inverse_transform](y_calc_sc)
	
	if roi_out == "no"
    	return y[:,1] - y_calc, y_calc
	elseif roi_out == "yes" 
		return y[:,1] - y_calc, y_calc, [interest_x interest_y]
	else
		error("roi_out should be set to yes or no")
	end
end
	