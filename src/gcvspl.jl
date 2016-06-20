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

function gcvspl(x::Array{Float64,1},y::Array{Float64,1},ese::Array{Float64,1},SmoothSpline::Float64;SplineOrder::Int32 = Int32(2),SplineMode::Int32 = Int32(3))
	"""
    
    c, wk, ier = gcvspline(x,y,ese,SplineSmooth,SplineOrder,SplineMode)

    INPUTS:
	
    x: Float64 Array, the independent variables
    y: Float64 Array, the observations (we assume here that you want to use this spline only on 1 dataset... see gcvspl.f if not)
    ese: Float64 Array, the errors on y
    SplineSmooth: Float64, the smoothing factor
    SplineOrder (M parameter in gcvspl.f): Int32, the half order of the required B-splines. default: splorder = 2 (cubic)
    SplineOrder = 1,2,3,4 correspond to linear, cubic, quintic, and heptic splines, respectively. 
    SplineMode (MD parameter in gcvspl.f) is the Optimization mode switch:
	    default:   SplineMode = 2 (General Cross Validated)
	               SplineMode = 1: Prior given value for p in VAL
	                         (VAL.ge.ZERO). This is the fastest
	                         use of GCVSPL, since no iteration
	                         is performed in p.
	               SplineMode = 2: Generalized cross validation.
	               SplineMode = 3: True predicted mean-squared error,
	                         with prior given variance in VAL.
	               SplineMode = 4: Prior given number of degrees of
	                         freedom in VAL (ZERO.le.VAL.le.N-M).
	               SplineMode  < 0: It is assumed that the contents of
	                         X, W, M, N, and WK have not been
	                         modified since the previous invoca-
	                         tion of GCVSPL. If MD < -1, WK(4)
	                         is used as an initial estimate for
	                         the smoothing parameter p.  At the
	                         first call to GCVSPL, MD must be > 0.
	               Other values for |MD|, and inappropriate values
	               for VAL will result in an error condition, or
	               cause a default value for VAL to be selected.
	               After return from MD.ne.1, the same number of
	               degrees of freedom can be obtained, for identical
	               weight factors and knot positions, by selecting
	               |MD|=1, and by copying the value of p from WK(4)
	               into VAL. In this way, no iterative optimization
	               is required when processing other data in Y. 

    OUPUTS:
	
		c: the spline coefficients
		WK: work vector, see gcvspl.f
		IER: error parameter. 
		IER = 0: Normal exit 
		IER = 1: M.le.0 .or. N.lt.2*M
		IER = 2: Knot sequence is not strictly
		increasing, or some weight
		factor is not positive.
		IER = 3: Wrong mode parameter or value.
		
    SEE gcvspl.f FOR MORE INFORMATION
	
	"""
	if size(y,2)>1
		error("This function does not accept multiple dataset. Please provide only one dataset at a time.")
		#elseif size(x,1) != size(y,1) || size(ese,1) != size(y,1) || size(x,1) != size(ese,1)
		#error("This function only accepts x, y and ese arrays of the same lengths.")
	end
	
	WX = 1. ./(ese.^2) # relative variance of observations
	WY = zeros([1])+1. # systematic errors... not used so put them to 1
	VAL = (SmoothSpline.*ese).^2

	# M = SplineOrder # No need to redefine but just this comment to keep record
	N = size(x,1) #same as length(y)
	K = 1 # number of y columns
	# MD = SplineMode # No need to redefine but just this comment to keep record
	NC = length(y)

	c = ones(N,NC)
	WK = ones(6*(N*SplineOrder+1)+N,1)
	IER=Int32[1]
	ccall( (:gcvspl_, abspath("../Dependencies/gcvspline/libgcvspl.so")), Void, (Ptr{Float64},Ptr{Float64},Ptr{Cint},Ptr{Float64},Ptr{Float64},Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Float64},Ref{Float64},Ptr{Cint},Ref{Float64},Ref{Cint}),x,y,&N,WX,WY,&SplineOrder,&N,&K,&SplineMode,VAL,c,&NC,WK,IER)
	return c, WK, IER
end

function splderivative(xfull::Array{Float64},xparse::Array{Float64},cparse::Array{Float64};SplineOrder::Int32 = Int32(2), L::Int32 = Int32(1), IDER::Int32 = Int32(0))
    """
    Wrapper to the SPLDER function of gcvspl.f, for interpolation purpose
    
	INPUTS:
	
	    xfull: Float64 Array, contains the entire x range where the spline has to be evaluated
	    xparse: Float64 Array, contains the x values of interpolation regions   
	    WARNING!!! => xparse[0] <= xfull[0] <= xparse[n] 
	    cparse: Float64 Array, is the evaluated spline coefficients returned by gcvspl for xparse
    
    OPTIONS:    
	
        splineorder (integer): is the spline order, 
        default: splineorder = 2 (cubic)
        L (integer): see gcvspl.f for details, default: L = 1
        IDER: the Derivative order required, with 0.le.IDER 
        and IDER.le.2*M. If IDER.eq.0, the function
        value is returned; otherwise, the IDER-th
        derivative of the spline is returned.
        
    SEE gcvspl.f FOR MORE INFORMATION
    """   
	# Note: In the following t, x and c in the fortran code are renamed  xfull, xparse and cparse for clarity
    N = Int32(size(xparse,1))    
    q = zeros(2*SplineOrder,1)  # working array
    y_calc::Array{Float64} = zeros(size(xfull,1),1) # Output array
    
    # we loop other xfull to create the output values
	for i =1:size(y_calc,1)
	    y_calc[i] = ccall( (:splder_, abspath("../Dependencies/gcvspline/libgcvspl.so")), Float64, (Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Cint},Ptr{Float64}),&IDER, &SplineOrder, &N, &xfull[i], xparse, cparse, &L, q)
	end
	return y_calc
end