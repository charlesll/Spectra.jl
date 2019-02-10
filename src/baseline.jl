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



Allows subtracting a baseline under a x y spectrum.

Parameters
----------
x_input : ndarray
	x values.
y_input : ndarray
	y values.
bir : ndarray
	Contain the regions of interest, organised per line. 
	For instance, roi = np.array([[100., 200.],[500.,600.]]) will 
	define roi between 100 and 200 as well as between 500 and 600.
	Note: This is NOT used by the "als" and "arPLS" algorithms, but still is a requirement when calling the function.
	bir and method probably will become args in a futur iteration of rampy to solve this.
methods : str
	"poly": polynomial fitting, with splinesmooth the degree of the polynomial.
	"unispline": spline with the UnivariateSpline function of Scipy, splinesmooth is 
				 the spline smoothing factor (assume equal weight in the present case);
	"gcvspline": spline with the gcvspl.f algorythm, really robust. 
				 Spectra must have x, y, ese in it, and splinesmooth is the smoothing factor;
				 For gcvspline, if ese are not provided we assume ese = sqrt(y). 
				 Requires the installation of gcvspline with a "pip install gcvspline" call prior to use;
	"exp": exponential background;
	"log": logarythmic background;
	"rubberband": rubberband baseline fitting;
	"als": automatic least square fitting following Eilers and Boelens 2005;
	"arPLS": automatic baseline fit using the algorithm from Baek et al. 2015 
			 Baseline correction using asymmetrically reweighted penalized least squares smoothing, Analyst 140: 250-257.

Options
------
polynomial_order : Int
	The degree of the polynomial (0 for a constant), default = 1.
s : Float
	spline smoothing coefficient for the unispline and gcvspline algorithms.
lam : Float
	float, the lambda smoothness parameter for the ALS and ArPLS algorithms. Typical values are between 10**2 to 10**9, default = 10**5.
p : Float
	float, for the ALS algorithm, advised value between 0.001 to 0.1, default = 0.01.
ratio : float
	ratio parameter of the arPLS algorithm. default = 0.01.
niter : Int
	number of iteration of the ALS algorithm, default = 10.
p0_exp : List
	containg the starting parameter for the exp baseline fit with curve_fit. Default = [1.,1.,1.].
p0_log : List
	containg the starting parameter for the log baseline fit with curve_fit. Default = [1.,1.,1.,1.].

Returns
-------
out1 : ndarray
	Contain the corrected signal.
out2 : ndarray
	Contain the baseline.

"""
function baseline(x::Array{Float64},y::Array{Float64},roi::Array{Float64},basetype::AbstractString;polynomial_order=1, s = 1.0, lam = 10^5, p = 0.01, ratio = 0.01, niter = 10, p0_exp = [1.,1.,1.],p0_log =[1.,1.,1.])
	yout, bas = rampy[:baseline](x,y,roi,basetype,polynomial_order=polynomial_order, s=s, lam=lam, niter=niter, p0_exp=p0_exp, p0_log=p0_log)
end
