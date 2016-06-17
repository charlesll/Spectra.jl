.. _Splines:
***********************
Splines
***********************

-------------------
Introductory notes
-------------------

Not all the splines packages provide the same performances for data smoothing and interpolation. By experience, the Dierckx spline package ("Dspline" option in the baseline() function) provides a good starting point, but is not as usefull as other spline packages.

The csaps function of Matlab uses the SMOOTH Fortran library, and provides better smoothing capabilities for noisy data. Similarly, the GCVSPL Fortran package from Woltring (1986) also provides a very robust way to smooth and interpolate noisy data.

This GCVSPL spline package is called directly by Julia (through a ccall()) in the baseline function, with the options of a cubic spline with least-square data fitting. The smoothing is done with scaling the variances of the data points (VAR variable in the GCVSPL.f package) that is provided to the GCVSPL.f program.

Now, while baseline() should be well suited for most users needs, it uses cubic splines that are not always the best answers to some problems. For instance, quadratic splines may be more robust in some cases. You can change that by providing the spline order to baseline() as SplOrder = 2 for instance.

In case you want to have even more control on GCVSPL.f, and use its internal tricks and tweeks, the following lines will provide you the documentation of the two functions allowing you to calculate the spline coefficients and to evaluate the spline values at specific x entries.

------------------------------
Function gcvspl
------------------------------

This function allows you to calculate the spline coefficients. It calls gcvspline subroutine of the program GCVSPL.f as::

    c, wk, ier = gcvspline(x,y,ese,SplineSmooth,SplineOrder,SplineMode)

INPUTS:
	
    x: Array{Float64}, the independent variables;
	
    y: Array{Float64}, the observations (we assume here that you want to use this spline only on 1 dataset... see gcvspl.f if not);
	
    ese: Array{Float64}, the errors on y;
	
    SplineSmooth: Float64, the smoothing factor;
	
    SplineOrder (M parameter in gcvspl.f): Int32, the half order of the required B-splines. default: splorder = 2 (cubic). SplineOrder = 1,2,3,4 correspond to linear, cubic, quintic, and heptic splines, respectively. 
	
    SplineMode (Int32, MD parameter in gcvspl.f) is the Optimization mode switch:
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
							 
	               Other values for SplineMode, and inappropriate values
	               for VAL will result in an error condition, or
	               cause a default value for VAL to be selected.
	               After return from MD.ne.1, the same number of
	               degrees of freedom can be obtained, for identical
	               weight factors and knot positions, by selecting
	               SplineMode=1, and by copying the value of p from WK(4)
	               into VAL. In this way, no iterative optimization
	               is required when processing other data in Y. 

OUPUTS:
	
	c: Array{Float64}, the spline coefficients;
	
	WK: Array{Float64}, working vector, see gcvspl.f;
	
	IER: error parameter
	 
		IER = 0: Normal exit 
		
		IER = 1: M.le.0 .or. N.lt.2*M
		
		IER = 2: Knot sequence is not strictly increasing, or some weight factor is not positive.
		
		IER = 3: Wrong mode parameter or value.
		
SEE GCVSPL.f and Woltring (1986) for even more information.

------------------------------
Function splderivative
------------------------------

After a call of gcvspl, this function allows you to calculate the spline values for given x entries. It is called as::

	splderivative(xfull,xparse,cparse;SplineOrder::Int32 = Int32(2), L::Int32 = Int32(1), IDER::Int32 = Int32(0))

INPUTS:

	xfull: Array{Float64}, contains the entire x range where the spline has to be evaluated;
	
	xparse: Array{Float64}, contains the x values of interpolation regions. WARNING!!! => xparse[0] <= xfull[0] <= xparse[n]; 
	
	cparse: Array{Float64}, is the evaluated spline coefficients returned by gcvspl for xparse.
    
OPTIONS:    
	
	splineorder (integer): is the spline order, 
		default: splineorder = 2 (cubic) (yes 2 is for a cubic spline, not 3. Be careful with this!)

		L (integer): see gcvspl.f for details, default: L = 1

		IDER: the Derivative order required, with 0.le.IDER and IDER.le.2*M. If IDER.eq.0, the function value is returned; otherwise, the IDER-th derivative of the spline is returned.
        
SEE GCVSPL.f and Woltring (1986) for even more information.