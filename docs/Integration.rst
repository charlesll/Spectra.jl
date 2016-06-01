**************
 Integration functions
**************

Functions are available in integration.jl to be able to integrate the area under a region of a spectrum, or to integrate the area under Gaussian, Lorentzian or other peaks built from the functions provided in functions.jl.

Available functions
At the moment, only two functions are available. trapz is an implementation in Julia of the trapezoïdal integration function available in Matlab for instance. Call it as:

    area = trapz(x, y)

x: vector containing the x values;
y: vector containing the y values;

This function is particularly helpful to calculate the area under a portion of a spectrum, and uses it for various purposes (normalisation, area comparison, etc.).

gaussianarea allows to calculate the area under a gaussian peak from its half-width at half maximum (hwhm) and its amplitude, with the possibility of calculating the error based on the inputs of the errors on hwhm and amplitude. Call it as:

    area, esearea = gaussianarea(Amplitude,HWHM::Array{Float64}; eseAmplitude::Array{Float64} = [0.0], eseHWHM::Array{Float64} = [0.0])

Amplitude: Array{Float64} containing the amplitudes
HWHM: Array{Float64} containing the half width at half maximum of the peaks
eseAmplitude: Array{Float64} containing the errors affecting Amplitude
eseHWHM: Array{Float64} containing the errors affecting HWHM

area = Array{Float64} contaning the areas
esearea = Array{Float64} contaning the propagated errors affecting the areas calculations. If eseAmplitude and eseHWHM are not indicated, the function will return 0.0 for esearea, the error made on the area calculation.

----------
 To Do
----------

gaussianarea will change in peakarea, with the option to choose between the shape of the peak (gaussian, lorentzian, etc.)

Should we think to add Simpson's rule integration also?
