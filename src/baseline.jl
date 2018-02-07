#############################################################################
#Copyright (c) 2016-2017 Charles Le Losq
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

"""
	baseline(x::Array{Float64},y::Array{Float64},roi::Array{Float64},basetype::AbstractString;p=1.0,lambda = 10^5,SplOrder=3,roi_out="no")

Baseline subtraction can be made with using the baseline function:

INPUTS:

	x: Array{Float64}

containing the x values;

	y: Array{Float64}

containing the y values;

	roi: Array{Float64}

containing the region of interest, i.e. the places where you want to fit the baseline. For instance, if the baseline should fit the regions comprised between 750 and 800 cm^{-1}, and 1250 and 1300 cm^{-1}: roi = [750. 800.; 1250. 1300.];

	basetype: AbstractString

the type of baseline that you want to use. For now, polynomial and cubic spline baselines are available. Indicate the type you want as:

Polynomial baseline: enter "poly" for basetype, then the polynomial degree as p.

Dierckx cubic spline baseline: enter "Dspline" for basetype, then the smoothing degree as p. This uses the Dierckx package: https://github.com/kbarbary/Dierckx.jl

Generalised Cross-Validated baseline: enter "gsvspline" for basetype, then the smoothing degree as p.

Kernel Ridge Regression: enter "KRregression" for basetype, no need to provide p.

Support Vector Machines regression: enter "SVMregression" for basetype, no need to provide p.

ALS algorithm: enter "als" for automatic baseline fitting following Eilers and Boelens (2005).

arPSL algorithm: enter "arPLS" for automatic baseline fitting following Baek et al. (2015).

whittaker algorithm: enter "whittaker" to use the whittaker smoother described in Eiler (2003) which fit the signal in the roi regions. See also function whitsmdd().

OPTIONS:

	p:: Float64

Default = 1.0. If using gcvspline or Dspline, this number indicates the spline smoothing coefficient. If using "poly", it is the degree of the polynomial function to be fitted. Please enter a float number (1.0, 2.0 or 3.0 for splines of order 1, 2 or 3), and it is automatically converted to an Integer for the polyfit function.

For the ALS algorithm, choose p in the range 0.001 - 0.1.

For the arPLS algorithm, p corresponds to the breaking ratio parameter in the paper of Baek et al. (2015). Test different values of p, starting at high p values (~0.1).

	lambda: Float64

smoothing parameter of the ALS, arPLS and whittaker algorithms, recommended values in the range 10^2 - 10^9; Default = 10^5.

	SplOrder: Integer

the spline coefficient to be used with the Dspline or gcvspline options. Default = 3.

	roi_out: String, "no" or "yes".

This will result in an additional output matrix containing the y signal in the roi regions of interest, which can then be used to plot and to evaluate the roi provided to the baseline function.

	niter: Int

number of iterations for the ALS algorithm. Default = 10.

OUTPUTS:

(are combined in a tuple in one array if only one output variable is provided)

	y_corr: Array{Float64}

the spectrum corrected from its baseline;

	bass: Array{Float64}

the baseline.

OPTIONAL OUTPUT:

	y_roi_out: Array{Float64}

an 2 column array containing the initial x-y pairs of the signal in the roi regions of interest.

NOTES:

Errors on measurenements are automatically provided as sqrt(y) in gcvspline. For further options, please use the gcvspl and splderivative functions that directly call the GCVSPL and SPLDER function of the gcvspl.f program (Holtring, 1986). Further informations for the use of splines are given in the Splines section, see :ref:`Splines`.

The Kernel Ridge and Support Vector Machines regression algorithms call the Scikit Learn library, available in Python. This library thus SHOULD be installed. They are machine learning algorithms that will try to automatically fit the baseline in the provided regions of interest. They are slower that splines, but have the advantage of avoiding the (sometimes painful) tuning of the spline coefficients.

The Kernel Ridge and Support Vector Machines regression algorithms used a Cross-Validated approach to increase the generalisation and avoid overfitting. The GridSearchCV function of SciKit Learn is called, with 5 fold cross-validation and the following gridsearch parameters:

- For KRregression:

	param_grid=Dict("alpha"=> [1e0, 0.1, 1e-2, 1e-3],"gamma"=> logspace(-4, 4, 9));

- For SVMregression:

	param_grid=Dict("C"=> [1e0, 1e1, 1e2, 1e3],"gamma"=> logspace(-4, 4, 9)).

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

    bas = baseline(x,y,roi,"Dspline",p=0.01)

p there is the smoothing parameter used.

"""

function baseline(x::Array{Float64},y::Array{Float64},roi::Array{Float64},basetype::AbstractString;p=1.0,lambda = 10^5,SplOrder=3,roi_out="no",niter = 10)
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
	# Getting mean and std for scaling
	x_mean = mean(interest_x)
	y_mean = mean(interest_y)

	x_std = std(interest_x)
	y_std = std(interest_y)

	# Scaling the data
	x_bas_sc = (interest_x-x_mean)/x_std
	y_bas_sc = (interest_y-y_mean)/y_std

	x_sc = (x-x_mean)/x_std
	y_sc = (y-y_mean)/y_std

	######## POLYNOMIAL BASELINE
	if basetype == "poly"
		degree = round(Int,p)
		if degree == 0
			y_calc_sc::Array{Float64} = ones(size(y,1),1).*mean(y_bas_sc)
		else
			best_p = polyfit(x_bas_sc[:,1],y_bas_sc[:,1],round(Int,p))
			y_calc_sc = polyval(best_p,x_sc[:,1])
		end
		y_calc = (y_calc_sc*y_std)+y_mean # unscalling

	######## DIERCKX SPLINE BASELINE
	elseif basetype == "Dspline"
		spl = Spline1D(x_bas_sc[:,1],y_bas_sc[:,1],s=p[1],bc="extrapolate",k=SplOrder)
		y_calc_sc = evaluate(spl,x_sc[:,1])
		y_calc = (y_calc_sc*y_std)+y_mean # unscalling

	######## GCV SPLINE BASELINE
	elseif basetype == "gcvspline"
		#c, WK, IER = gcvspl_julia(x_bas_sc[:,1],y_bas_sc[:,1],ese_interest_y[:,1],p[1];SplineOrder = Int32(SplOrder-1)) # with cubic spline as the default
		#y_calc_sc = splderivative_julia(x_sc[:,1],x_bas_sc[:,1],c,SplineOrder= Int32(SplOrder-1))

		# implementation using the Python wrapper of gcvspline.
		ese_interest_y = ones(size(y_bas_sc))
		#flt = pygcvspl[:MSESmoothedNSpline](x_bas_sc[:,1],y_bas_sc[:,1],1./(p[1].*ese_interest_y[:,1])) #variance_metric=p[1].^2
		flt = pygcvspl[:SmoothedNSpline](x_bas_sc[:,1],y_bas_sc[:,1],p[1])
		y_calc_sc = flt(x_sc)
		y_calc = (y_calc_sc*y_std)+y_mean # unscalling

	######## WHITTAKER BASELINE (DIRECT RETURN)
	elseif basetype == "whittaker"

		w = zeros(size(y,1))
		w[interest_index] = 1.0

		y_calc = whitsmdd(x,y,w,lambda)


	######## ALS AUTOMATIC BASELINE
	elseif basetype == "als"
		# als algorithm from Matlab code in Eilers et Boelens 2005

		# Estimate baseline with asymmetric least squares
		m = length(y)
		D = diff(diff(speye(m))) # use if data equally spaced
		#D = ddmat(x, 2) # from Eilers 2003: use ddmat if data are not equally spaced...
		#w = ones(N, 1)
		w = zeros(size(y,1)) # Modification following email of Baek, thanks!
		w[interest_index] = 1.0

		for it = 1:niter
			W = spdiagm((w,), 0, m, m)
			#C = chol(W + lambda * D' * D) # for whatever reason cholfac does not work well... so we solve the problem a bit differently from the matlab code. Results seem similar.
			global z = (W + lambda * D' * D) \ (w .* y) #C \ (C' \ (w .* y))
			w[w.!=0] = p * (y[w.!=0] .> z[w.!=0]) + (1 - p) * (y[w.!=0] .< z[w.!=0])
		end

		y_calc = collect(z)

	######## arPLS AUTOMATIC BASELINE
	elseif basetype == "arPLS"
		# arPLS algorithm of Baek et al. 2015 Analyst 140: 250-257

		#w = ones(N, 1)
		w = zeros(size(y,1)) # Modification following email of Baek, thanks!
		w[interest_index] = 1.0

		# Estimate baseline with asymmetric least squares
		N = length(y)
		#D = ddmat(x, 2) # from Eilers 2003: use ddmat if data are not equally spaced...
		D = diff(diff(speye(N))) # use if data equally spaced
		H = lambda * D' * D

		while true
			W = spdiagm((w,), 0, N, N)
			global z = (W + H) \ (w .* y) # We don't use the exact same Cholesky decomposition as in the code but the \ operator
			d = y - z
			# make d- and get w^t with m and s
			dn = d[d.<0]
			m = mean(dn)
			s = std(dn)
			wt = 1./(1 + exp( 2* (d-(2*s-m))/s ) )
			# check exit condition and backup
			if norm(w-wt)/norm(w) .< p; break; end
			#w = wt;
			w[w.!=0]=wt[w.!=0] # Modification following email of Baek, thanks!
		end
		y_calc = collect(z)

	######## KERNEL RIDGE REGRESSION WITH SCIKIT LEARN
	elseif basetype == "KRregression"

		clf = kernel_ridge[:KernelRidge](kernel="rbf", gamma=0.1)
		kr = model_selection[:GridSearchCV](clf,cv=5,param_grid=Dict("alpha"=> [1e1, 1e0, 0.5, 0.1, 5e-2, 1e-2, 5e-3, 1e-3,1e-4],"gamma"=> logspace(-4, 4, 9)))# GridSearchCV for best parameters
		kr[:fit](reshape(x_bas_sc,size(x_bas_sc,1),1), squeeze(reshape(y_bas_sc[:,1],size(y_bas_sc,1),1),2)) #SciKit learn is expecting a y vector, not an array...
		y_calc_sc = kr[:predict](reshape(x_sc,size(x_sc,1),1))
		y_calc = (y_calc_sc*y_std)+y_mean # unscalling

	######## SUPPORT VECTOR MACHINES REGRESSION WITH SCIKIT LEARN
	elseif basetype == "SVMregression"

		clf = svm[:SVR](kernel="rbf", gamma=0.1)
		svr = model_selection[:GridSearchCV](clf,cv=5,param_grid=Dict("C"=> [1e-1, 1e0, 1e1, 1e2, 1e3],"gamma"=> logspace(-4, 4, 9))) # GridSearchCV for best parameters
		svr[:fit](reshape(x_bas_sc,size(x_bas_sc,1),1), squeeze(reshape(y_bas_sc[:,1],size(y_bas_sc,1),1),2)) #SciKit learn is expecting a y vector, not an array...
		y_calc_sc = svr[:predict](reshape(x_sc,size(x_sc,1),1))
		y_calc = (y_calc_sc*y_std)+y_mean # unscalling

	######## RAISING ERROR IF NOT THE GOOD CHOICE
	else
        error("Not implemented, choose between poly, Dspline, gcvspline, als, arPLS, KRregression, SVMregression.") # for now GPregression is hidden
    end

	if roi_out == "no"
    	return y[:,1] - y_calc, y_calc
	elseif roi_out == "yes"
		return y[:,1] - y_calc, y_calc, [interest_x interest_y]
	else
		error("roi_out should be set to yes or no")
	end
end

"""
    ddmat(x, d)

Compute divided differencing matrix of order d

Input:

    x:  vector of sampling positions

    d:  order of diffferences
Output

    D:  the matrix; D * Y gives divided differences of order d

Matlab version in Paul Eilers, 2003

Julia translation by Charles Le Losq 2017
"""

function ddmat(x, d)


    m = length(x)
    if d == 0
        D = speye(m)
    else
        dx = x[(d + 1):m] - x[1:(m - d)]
        V = spdiagm((1 ./ dx,), 0, m - d, m - d)
        D = V * diff(ddmat(x, d - 1))
    end

    return D
end



"""
    whitsmdd(x,y,w,lambda;d=2)
Whittaker smoother with divided differences (arbitrary spacing of x)

Input:

  x:      data series of sampling positions (must be increasing)

  y:      data series, assumed to be sampled at equal intervals

  lambda: smoothing parameter; large lambda gives smoother result

  d:      order of differences (default = 2)

Output:

  z:      smoothed series

Matlab version by Paul Eilers, 2003
Julia translation by Charles Le Losq 2017

"""
function whitsmdd(x,y,w,lambda;d=2)

    # Smoothing
    m = length(y)
    E = speye(m)
    D = ddmat(x, d)
    W = spdiagm((w,), 0, m, m);
    # avoiding the col implementation in julia, direct use of \
    #C = chol(W + lambda * D' * D);
    #z = C \ (C' \ (w .*y));
    z = (W + lambda * D' * D) \ (w .* y)
    return z
end
