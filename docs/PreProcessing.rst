**************
 Pre-Processing Spectra
**************

----------
"Long" correction for glass Raman data
----------

The correction of glasses Raman spectra from temperature and excitation line effects (see REF TO ADD) can be done using the function:

    x, long, eselong = long(data::Array{Float64},temp::Float64,wave::Float64)

INPUTS:

data:
temp:
wave:

OUTPUTS:

(are combined in one array if only one output name is given)

x: array{Float64} containing the x values;
long: array{Float64} containing the corrected y values;
eselong: array{Float64} containing the errors calculated as sqrt(y) on raw data and propagated after the correction.

----------
Baseline subtraction function
----------

Baseline subtraction can be made with using the baseline function:

    y_corr, bass = baseline(x::Array{Float64},y::Array{Float64},roi::Array{Float64},basetype::AbstractString,p::Array{Float64})

INPUTS:

x: array{Float64} containing the x values;
y: array{Float64} containing the y values;
roi: an array containing the "region of interest", i.e. the places where you want to fit the baseline. For instance, if the baseline should fit the regions comprised between 750 and 800 cm^{-1}, and 1250 and 1300 cm^{-1}:

    roi = [[750. 800.]; [1250. 1300.]]

basetype: the type of baseline that you want to use. For now, polynomial and cubic spline baselines are available. Indicate the type you want as:

Polynomial baseline: enter "poly" for basetype, then the polynomial degree as p.

Dierckx cubic spline baseline: enter "Dspline" for basetype, then the smoothing degree as p.

Generalised Cross-Validated baseline: enter "gsvspline" for basetype, then the smoothing degree as p. 
Errors on measurenements are automatically provided as sqrt(y) in gcvspline. For further options, please use the gcvspl and splderivative functions that directly call the GCVSPL and SPLDER function of the gcvspl.f program (Holtring, 1986).

OUTPUTS:

(are combined in one array if only one output variable is provided)

y_corr: array{Float64}, the spectrum corrected from its baseline;
bass: array{Float64}, the baseline.

----------
Examples
----------

For instance, for subtracting a constant baseline between 1250 and 1300 cm^{-1}:

    roi = [[1250. 1300.]]
    basetype = "poly"
    p = [1.0]

such that:

    p = [1.0]
    bas = baseline(x,y,roi,"poly",p)

and the corrected spectrum is

    y_corr = y - bas

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

----------
Tips
----------

Always be careful to enter float numbers! The function will return an error if you do not do that.

For the spline, do not hesitate to test a broad range in term of order of magnitudes for the smoothing parameter!

----------
To Do
----------

- For now, x and y should contain only one column (one dataset at a time). In the futur, an option allowing to fit entire dataset may be provided?
- Adding access to the SMOOTH spline library, used by the csaps Matlab function.

----------
References
----------

Woltring, 1986, A FORTRAN package for generalized, cross-validatory spline smoothing and differentiation. Adv. Eng. Softw. 8:104-113. 
