**************
 Baseline subtraction
**************

----------
Baseline subtraction function
----------

Baseline subtraction can be made with using the baseline function:

    y_corr, bass = baseline(x::Array{Float64},y::Array{Float64},roi::Array{Float64},basetype::AbstractString,p::Array{Float64})

x: array{Float64} containing the x values;
y: array{Float64} containing the y values;
roi: an array containing the "region of interest", i.e. the places where you want to fit the baseline. For instance, if the baseline should fit the regions comprised between 750 and 800 cm^{-1}, and 1250 and 1300 cm^{-1}:

    roi = [[750. 800.]; [1250. 1300.]]

basetype: the type of baseline that you want to use. For now, polynomial and cubic spline baselines are available. Indicate the type you want as:

Polynomial baseline: enter "poly" for basetype, then the polynomial degree as p.
Cubic spline baseline: enter "spline" for basetype, then the smoothing degree as p.

outputs: y_corr and bass are respectively the spectrum corrected from its baseline, and the baseline.

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
    basetype = "spline"
    s = [0.01]
    bas = baseline(x,y,roi,basetype,s)

s there is the smoothing parameter used. The cubic spline uses the Dierckx package initially written in Fortran and used in Julia: https://github.com/kbarbary/Dierckx.jl

It gives pretty good results, but in the future I will wrap the GCVSPL spline, that may be more robust for very noisy data from my experience under using both algorithms under Python. GCVSPL is already available in RamPy if you really needs it. Another solution would be to wrap the csaps function of Matlab that gives pretty good results too on noisy data.

----------
Tips
----------

Always be careful to enter float numbers! The function will return an error if you do not do that.

For the spline, do not hesitate to test a broad range in term of order of magnitudes for the smoothing parameter!

----------
To Do
----------
For now, x and y should contain only one column (one dataset at a time). In the futur, an option allowing to fit entire dataset will be provided.
GCV splines are also going to be added. The Fortran code is already provided with Spectra.jl. I just need to wrap it to provide this functionality.
