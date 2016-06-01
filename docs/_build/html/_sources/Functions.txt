**************
 Available Functions
**************

The module functions.jl aims to contain several functions that are frequently used when dealing with spectra. For instance, it actually contains function for generating Gaussian peaks, polynomial curves, or lorentzian peaks.

----------
 Function poly
----------

This function just allows to build a polynomial curve as

    y = poly(p,x)

p: Vector{Float64} containing the polynomial parameters. For a linear curve, p = [1.0,1.0], for a second order polynomial, p = [1.0,1.0,1.0], etc.;
x: Array{Float64} containing the x values for calculation.

Output:

y: Array{Float64} containing the result of calculation.

----------
 Function gaussiennes
----------

gaussiennes, written in the plural french form there, is a function that allows to build gaussian peaks.

The gaussian function used there is:

    y = amplitude x exp(-ln(2) x [(x-centre)/hwhm]^2 )

To call gaussienne:

    gaussiennes(amplitude,centre,hwhm,x;style)

Inputs:
amplitude: Array{Float64} containing the peaks amplitudes;
centre: Array{Float64} containing the peaks centres;
hwhm: Array{Float64} containing the peaks half-width at middle heights (hwhm);
x: Array{Float64} containing the x axis values;

Optional arguments:
style: ASCIIString = "None", see examples below.

Outputs:
y: Array{Float64} containing the calculated y values.

You can enter the amplitude, centre and half-width at half-maximum (hwhm) values as arrays of float 64 (even containing one float value), without specifying style. hwhm is proportional to the standard deviation sigma:

hwhm= sqrt(2xln(2)) x sigma
that is used in a normal distribution (see function normal_dist)

----------
 Examples
----------

To have four gaussian peaks centered at 800, 900, 1000 and 1100 cm-1 with hwhm of 50 cm-1 on a Raman spectrum, you will enter:

   y_calc, y_peaks = gaussiennes([1.0,1.0,1.0,1.0], [800.0,900.0,1000.0,1100.0], [50.0,50.0,50.0,50.0], x)

and y_peaks will contain in 4 columns the 4 different y values of the peaks, and y_calc their sum (the total model). Now, if you want to calculate more complex models, such as for instance contructing how the Raman peaks of water vary with pressure, you might like to parametrize the variations of the peak parameters rather than just fitting each spectrum. This will provide more robust fits of the spectra, as you will fit them together, and will also force you to find the correct underlying mathematical assumption.

The gaussiennes function allows you to do that. If you specify style = "poly", you can enter arrays for the amplitudes, centres and half-widths at half-maximum (hwhm) of the peaks, with in each column the coefficients for the polynomial variations of this parameters. The second column of x will need to contain the second variable for those polynomial functions.

Let's say for instance that we have one peak at 900 cm-1 in a pure material. It's frequency seems to linearly shift with increasing the amount of hydrogen in this material, but it's intensity is non-linearly increasing, following a quadratic variation. It's width seems constant.

How to write that with gaussiennes? Well, first you need to construct a relevant x axis: first column contains the frequency, and the second one contains the chemical variable value. In our case, we want to model the peak between 800 and 1000 cm-1, for 1 wt% H. So we have an x array build like

    frequency = collect(800:1:1000)
    x = ones(length(frequency),2)
    x[:,1] = frequency[:]
    x[:,2] = 1.0

Ok, now lets build our y peaks:

    amplitudes = [1.0 0.1 0.1]
    frequencies = [900.0 2.0]
    hwhm = 20.0

    y_calc, y_peaks = gaussiennes(amplitudes, frequencies, hwhm, x)

This should provide you how the shape of the peak is as a function of both the frequency and the chemical composition there. If you want to go further, you might just want to stick gaussiennes in a loop, and play with creating various peaks with changing the chemical parameter in the x[:,2] column!

----------
 Functions lorentziennes, pseudovoigts, pearson7, normal_dist:
----------

Those functions return peaks that follow lorentzian, Pseudo-voigts, Pearson7 distributions. normal_dist return a peak following the normal distribution. They are implemented exactly the same way as gaussiennes, except that they return different peak shapes. We therefore refer the reader to the doc of the gaussiennes function for their uses. Inputs for those functions are:

    lorentziennes(amplitude::Array{Float64},centre::Array{Float64},hwhm::Array{Float64},x::Array{Float64};style::ASCIIString = "None")
    pearson7(a1::Array{Float64},a2::Array{Float64},a3::Array{Float64},a4::Array{Float64},x::Array{Float64};style::ASCIIString = "None")
    pseudovoigts(amplitude::Array{Float64},centre::Array{Float64},hwhm::Array{Float64},ps_fraction::Array{Float64},x::Array{Float64};style::ASCIIString = "None")
    normal_dist(nd_amplitudes::Array{Float64},nd_centres::Array{Float64},nd_sigmas::Array{Float64},x::Array{Float64})
