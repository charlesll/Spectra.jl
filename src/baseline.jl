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

include("gcvspl_wrapper.jl")

function baseline(x::Array{Float64},y::Array{Float64},roi::Array{Float64},basetype::AbstractString;p=1.0,SplOrder=3)
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
	ese_interest_y = sqrt(abs(interest_y)) # we assume errors as sqrt(y)

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
		best_p::Array{Float64} = polyfit(interest_x,interest_y,round(Int,p))
        y_calc::Array{Float64} = poly(best_p,x)
        return y[:,1] - y_calc, y_calc

	######## DIERCKX SPLINE BASELINE
	elseif basetype == "Dspline"
		spl = Spline1D(interest_x,interest_y,s=p[1],bc="extrapolate",k=SplOrder)
		y_calc = evaluate(spl,x[:,1])
		return y[:,1] - y_calc, y_calc

	######## GCV SPLINE BASELINE
	elseif basetype == "gcvspline"
		c, WK, IER = gcvspl_julia(interest_x[:,1],interest_y[:,1],ese_interest_y[:,1],p[1];SplineOrder = Int32(SplOrder-1)) # with cubic spline as the default
		y_calc = splderivative_julia(x[:,1],interest_x,c,SplineOrder= Int32(SplOrder-1))
		return y[:,1] - y_calc, y_calc

	######## KERNEL RIDGE REGRESSION WITH SCIKIT LEARN
	elseif basetype == "KRregression"
		# to be sure everything is in the right shape for SciKit Learn API
		interest_x = reshape(interest_x,size(interest_x,1),1)
		interest_y = reshape(interest_y,size(interest_y,1),1)
		
		# initialising the preprocessor scaler
		X_scaler = preprocessing[:StandardScaler]()
		Y_scaler = preprocessing[:StandardScaler]()

		X_scaler[:fit](interest_x)
		Y_scaler[:fit](interest_y)

		#scaling the data
		x_bas_sc = X_scaler[:transform](interest_x)
		y_bas_sc = Y_scaler[:transform](interest_y)
		x_sc = X_scaler[:transform](reshape(x,size(x,1),1))
		
		# constructing a GridSearchCV instance for grabing the best parameters
		clf = kernel_ridge[:KernelRidge](kernel="rbf", gamma=0.1)
		kr = grid_search[:GridSearchCV](clf,cv=5,param_grid=Dict("alpha"=> [1e1, 1e0, 0.5, 0.1, 5e-2, 1e-2, 5e-3, 1e-3,1e-4],"gamma"=> logspace(-4, 4, 9)))
		kr[:fit](x_bas_sc, squeeze(y_bas_sc,2)) #SciKit learn is expecting a y vector, not an array...
		y_kr_sc = kr[:predict](x_sc)
		y_kr = Y_scaler[:inverse_transform](y_kr_sc)
		return y[:,1] - y_kr, y_kr
		
	######## SUPPORT VECTOR MACHINES REGRESSION WITH SCIKIT LEARN
	elseif basetype == "SVMregression"
		# to be sure everything is in the right shape for SciKit Learn API
		interest_x = reshape(interest_x,size(interest_x,1),1)
		interest_y = reshape(interest_y,size(interest_y,1),1)
	
		# initialising the preprocessor scaler
		X_scaler = preprocessing[:StandardScaler]()
		Y_scaler = preprocessing[:StandardScaler]()

		X_scaler[:fit](interest_x)
		Y_scaler[:fit](interest_y)

		#scaling the data
		x_bas_sc = X_scaler[:transform](interest_x)
		y_bas_sc = Y_scaler[:transform](interest_y)
		x_sc = X_scaler[:transform](reshape(x,size(x,1),1))
	
		# constructing a GridSearchCV instance for grabing the best parameters
		clf = svm[:SVR](kernel="rbf", gamma=0.1)
		svr = grid_search[:GridSearchCV](clf,cv=5,param_grid=Dict("C"=> [1e-1, 1e0, 1e1, 1e2, 1e3],"gamma"=> logspace(-4, 4, 9)))
		svr[:fit](x_bas_sc, squeeze(y_bas_sc,2)) #SciKit learn is expecting a y vector, not an array...
		y_svr_sc = svr[:predict](x_sc)
		y_svr = Y_scaler[:inverse_transform](y_svr_sc)
		return y[:,1] - y_svr, y_svr
		
		######## GRADIENT BOOSTING ENSEMBLE REGRESSION WITH SCIKIT LEARN
		elseif basetype == "GPregression"
			# to be sure everything is in the right shape for SciKit Learn API
			interest_x = reshape(interest_x,size(interest_x,1),1)
			interest_y = reshape(interest_y,size(interest_y,1),1)
			ese_interest_y = reshape(ese_interest_y,size(ese_interest_y,1),1)
	
			# initialising the preprocessor scaler
			X_scaler = preprocessing[:StandardScaler]()
			Y_scaler = preprocessing[:StandardScaler]()

			X_scaler[:fit](interest_x)
			Y_scaler[:fit](interest_y)

			#scaling the data
			x_bas_sc = X_scaler[:transform](interest_x)
			y_bas_sc = Y_scaler[:transform](interest_y)
			x_sc = X_scaler[:transform](reshape(x,size(x,1),1))
	
			# constructing a GridSearchCV instance for grabing the best parameters
			gpr = gaussian_process[:GaussianProcess](corr="squared_exponential", theta0=1e-1,thetaL=1e-3, thetaU=1,nugget=(ese_interest_y[:] ./ interest_y[:]).^2,random_start=100)
			#svr = grid_search[:GridSearchCV](clf,cv=5,param_grid=Dict("C"=> [1e-1, 1e0, 1e1, 1e2, 1e3],"gamma"=> logspace(-2, 2, 5)))
			gpr[:fit](interest_x, squeeze(interest_y,2)) #SciKit learn is expecting a y vector, not an array...
			y_gpr = gpr[:predict](x_sc)
			#y_gpr = Y_scaler[:inverse_transform](y_gpr_sc)
			return y[:,1] - y_gpr, y_gpr
		
	######## RAISING ERROR IF NOT THE GOOD CHOICE
	else
        error("Not implemented, choose between poly, Dspline, gcvspline and KRregression.")
    end
end
