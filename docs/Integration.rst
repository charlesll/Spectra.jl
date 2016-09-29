**********************
 Integration functions
**********************

Spectra.jl provides functions that allow one to integrate the area under a region of a spectrum, or to calculate the area under Gaussian, Lorentzian or other bands.

---------------
 Function trapz
---------------

At the moment, only two functions are available. trapz is an implementation in Julia of the trapezoïdal integration function available in Matlab for instance. Call it as:

    area = trapz(x, y)

INPUTS:

	x: Vector{Float64} containing the x values;
	
	y: Vector{Float64} containing the y values;

OUTPUTS: 

	area: Float64, the trapezoidal integration value

This function is particularly helpful to calculate the area under a portion of a spectrum, and can be used for various purposes (normalisation, area comparison, etc.).

----------------------
 Function bandarea
----------------------

This function replaces the function gaussianarea in the version <0.1.9 of Spectra.jl. It allows to calculate the area under a specific band, with different shapes. For now, only Gaussian bands are supported, but other band shapes will be added soon. (This explains why gaussianarea is deprecated in favor of a more generic function)

gaussianarea allows to calculate the area under a gaussian peak from its half-width at half maximum (hwhm) and its amplitude, with the possibility of calculating the error based on the inputs of the errors on hwhm and amplitude. Call it as:

    area, esearea = band(Amplitude,HWHM; peak_shape = "Gaussian", error_switch = "no", eseAmplitude = [0.0], eseHWHM = [0.0])

INPUTS:

	Amplitude: Array{Float64}, contains the amplitudes (intensity) of the band(s);

	HWHM: Array{Float64}, contains the half width at half maximum of the peaks;

OPTIONS

	peak_shape: String, indicates the shape of the component. Only "Gaussian" is supported for now;
	
	error_switch: String, should be "yes" or "no". If "yes", the arrays containing the errors affecting the band amplitude and widhts should be provided in eseAmplitude and eseHWHM (see below);
	
	eseAmplitude: Array{Float64}, an array that contains the errors affecting Amplitude;

	eseHWHM: Array{Float64}, an array that contains the errors affecting HWHM;

OUTPUTS: 

	area: Array{Float64}, an array that contains the areas;

	if error_switch is set to "yes", then a second output is provided:
	
	esearea: Array{Float64}, an array that contains the propagated errors affecting the areas calculations.
