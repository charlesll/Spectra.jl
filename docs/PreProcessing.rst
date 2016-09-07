***********************
 Pre-Processing Spectra
***********************

------------------------------------------------------------
Temperature and frequency corrections for Raman spectra
------------------------------------------------------------

Raman spectra can be corrected from temperature and excitation line effects using the function

    x, y_corrected, ese_corrected = tlcorrection(data,temp,wave;correction="long",normalisation="area",density=2210.0)

INPUTS:

	data: Array{Float64}, input spectrum with x and y in first and second columns respectively;

	temp: Float64, the temperature in °C;

	wave: Float64, the wavenumber at which the spectrum was acquirred in nm.

OPTIONS:

	correction: String, the equation used for the correction. Choose between "long", "galeener", or "hehlen". Default = "long".
	
	normalisation: String, indicate if you want to normalise your signal or not. Choose between "intensity", "area", or "no". Default = "area".
	
	density: Float64, the density of the studied material in kg m-3, to be used with the "hehlen" equation. Default = 2210.0 (density of silica).

OUTPUTS:

(are combined in one array if only one output name is given)

	x: Array{Float64}, containing the x values;

	long: Array{Float64}, containing the corrected y values;

	eselong: Array{Float64}, containing the errors calculated as sqrt(y) on raw data and propagated after the correction.
	
NOTES:

This correction uses the formula reported in Galeener and Sen (1978), Mysen et al. (1982), Brooker et al. (1988) and Hehlen et al. (2010).

The "galeener" equation is the exact one reported in Galeener and Sen (1978), which is a modification from Shuker and Gammon (1970) for accounting of (vo - v)^4 dependence of the Raman intensity. See also Brooker et al. (1988) for further discussion.

The "long" equation is that of Galeener and Sen (1978) corrected by a vo^3 coefficient for removing the cubic meter dimension of the equation of "galeener". This equation has been used in Mysen et al. (1982), Neuville and Mysen (1996) and Le Losq et al. (2012).

The "hehlen" equation is that reported in Hehlen et al. (2010). It actually originates before this publication (Brooker et al. 1988). It uses a different correction that avoid crushing the signal below 500 cm-1. THerefore, it has the advantage of keeping intact the Boson peak signal in glasses.

------------------------------
Baseline subtraction function
------------------------------

Baseline subtraction can be made with using the baseline function:

    y_corr, bass = baseline(x,y,roi,basetype;p=1.0,SplOrder=3)

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
	
OUTPUTS:

(are combined in a tuple in one array if only one output variable is provided)

	y_corr: Array{Float64}, the spectrum corrected from its baseline;
	
	bass: Array{Float64}, the baseline.

NOTES:

Errors on measurenements are automatically provided as sqrt(y) in gcvspline. For further options, please use the gcvspl and splderivative functions that directly call the GCVSPL and SPLDER function of the gcvspl.f program (Holtring, 1986). Further informations for the use of splines are given in the Splines section, see :ref:`Splines`.

The Kernel Ridge and Support Vector Machines regression algorithms call the Scikit Learn library, available in Python. This library thus SHOULD be installed. They are machine learning algorithms that will try to automatically fit the baseline in the provided regions of interest. They are slower that splines, but have the advantage of avoiding the (sometimes painful) tuning of the spline coefficients.

The Kernel Ridge and Support Vector Machines regression algorithms used a Cross-Validated approach to increase the generalisation and avoid overfitting. The GridSearchCV function of SciKit Learn is called, with 5 fold cross-validation and the following gridsearch parameters:

	- For KRregression: param_grid=Dict("alpha"=> [1e0, 0.1, 1e-2, 1e-3],"gamma"=> logspace(-4, 4, 9));
	- For SVMregression: param_grid=Dict("C"=> [1e0, 1e1, 1e2, 1e3],"gamma"=> logspace(-4, 4, 9)).
	
Please see the SciKit Learn documentation at http://scikit-learn.org/stable/index.html for further details on the implementation of those technics, together with the source code of Spectra.jl.

----------
Examples
----------

For instance, for subtracting a constant baseline between 1250 and 1300 cm^{-1}:

    roi = [1250. 1300.]
	
    basetype = "poly"
	
    p = [1.0]
	
    y_corr, bas = baseline(x,y,roi,"poly",p)
	

For a linear baseline,

    p = [1.0, 1.0]
	
    bas = baseline(x,y,roi,"poly",p)

For a second order polynomial baseline,

    p = [1.0, 1.0, 1.0]
	
    bas = baseline(x,y,roi,"poly",p)

with the last coefficient will be the one in front of x^2. This can continue as you want by adding more 1.0 values to p.

For a cubic spline baseline fitting the basis of a peak centered at 1100 cm$^{-1}$ and with basis at 900 and 1250 cm^{-1}:

    roi = [[890. 910.]; [1250. 1300.]]
	
    basetype = "Dspline"
	
    s = [0.01]
	
    bas = baseline(x,y,roi,basetype,s)

s there is the smoothing parameter used. The cubic spline uses the Dierckx package initially written in Fortran and used in Julia: https://github.com/kbarbary/Dierckx.jl

---------------------------
Frequency shifts correction
---------------------------

In case your spectra are shifted from a reference value, Spectra offers several functions that allows you to correct it from this shift.

To correct a spectrum from a shift of P wavenumbers, you can simply call:

	xshift_direct(original_x::Array{Float64}, original_y::Array{Float64}, p::Float64)

INPUTS:

	full_x: Array{Float64}, the entire X axis of your spectrum
	
	original_x: Array{Float64}, the (shifted) X axis of original_y
	
	original_y: Array{Float64}, the signal
	
	p: Float64, the value of the shift
	
OUTPUTS:

	original_x: Array{Float64}, the entire X axis of your spectrum, we output it for records
	
	corrected_y: Array{Float64}, the signal corrected from the X shift.
	
	
This function uses the Dierckx spline to interpolate your signal after the correction of the shift.

Sometime, two signals from the same mineral show a shift in the X axis, while they share a common X axis. To correct from such thing, you can use the function:

	xshift_correction(full_x, full_shifted_y, ref_x, ref_y, shifted_y)
	
INPUTS:

	full_x: Array{Float64}, the entire X axis of your spectrum
	
	full_shifted_y: Array{Float64}, the entire shifted signal
	
	ref_x: Array{Float64}, the X axis of the reference ref_y signal
	
	ref_y: Array{Float64}, the reference signal
	
	shifted_y: Array{Float64}, the shifted signal
	
OUTPUTS:

	full_x: Array{Float64}, the entire X axis of your spectrum, we output it for records
	
	corrected_y: Array{Float64}, the signal corrected from the X shift
	
ref_x is the common X axis of two particular ref_y and shifted_y signals, that should be for instance an intense and well defined peak in your spectra. If ref_y and shifted_y do not share the same X axis, you can use first the Dierckx spline to re-sample one of them and have both sharing a common X axis. See the examples for further details.

--------------
References
--------------

Shuker, Reuben, and Robert Gammon. 1970. “Raman-Scattering Selection-Rule Breaking and the Density of States in Amorphous Materials.” Physical Review Letters 25 (4): 222–25.

Galeener, F. L., and Sen, P. N. 1978. “Theory of the First-Order Vibrational Spectra of Disordered Solids.” Physical Review B 17 (4): 1928–33.

Mysen, B. O., L. W. Finger, D. Virgo, and F. A. Seifert. 1982. “Curve-Fitting of Raman Spectra of Silicate Glasses.” American Mineralogist 67: 686–95.

Brooker et al. 1988 Assessment of correction procedures for reduction of Raman spectra. Journal of Raman Spectroscopy 19(2), 71-78.

Neuville, D. R., and B. O. Mysen. 1996. “Role of Aluminium in the Silicate Network: In Situ, High-Temperature Study of Glasses and Melts on the Join SiO₂-NaAl0₂.” Geochimica et Cosmochimica Acta 60: 1727–37.

Le Losq, C., D. R. Neuville, R. Moretti, and J. Roux. 2012. “Determination of Water Content in Silicate Glasses Using Raman Spectrometry: Implications for the Study of Explosive Volcanism.” American Mineralogist 97 (5-6): 779–90. doi:10.2138/am.2012.3831.

Hehlen, B. 2010. “Inter-Tetrahedra Bond Angle of Permanently Densified Silicas Extracted from Their Raman Spectra.” Journal of Physics: Condensed Matter 22 (2): 025401.
