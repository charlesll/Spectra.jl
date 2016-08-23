**************
 Peak fitting
**************

----------------
Model adjustment
----------------

Peak fitting is done with the JuMP framework (https://jump.readthedocs.org/en/latest/). Spectra.jl actually does not provide any peak fitting capacities by itself, but the combination of its functionality with JuMP helps making fitting procedure quite easy. An example is visible in the example section of Spectra.jl. 

One goal of Spectra is to promote the use of global optimisation models, where peak parameters are actually calculated from variation in other parameters (chemistry, temperature, etc.), or are shared between several spectra. I will provide very soon an example of such an approach. It be implemented in a few lines of code with combining Spectra and JuMP, and has the advantage of greatly reducing the errors of the fits.

------------------------------------
error calculation with bootstrapping
------------------------------------

Error calculation can be done with using bootstrapping. Spectra provides a function that allows generating K new datasetes, by resampling the existing dataset in a non-parametric or parametric way. 

The bootstrap function
----------------------

The bootstrap data generation function is called as:

    b_x_f, b_y_f = bootstrap(x::Array{Float64}, y::Array{Float64};boottype = "np")
	
INPUTS

	x: Array{Float64}, the x axis. It can have multiple columns.

	y: Array{Float64}, the y axis. It can have multiple columns.
	
OPTIONS

	boottype: ASCIIString, either "np" or "p", this is the type of bootstrapping performed. "np" performes a non-parametric resampling fo the dataset with replacement. "p" performs a parametric resampling. The data are resample from a gaussian distribution centered on the y values with errors that should be provided in the ese variable.
	
	ese: Array{Float64}, containing the errors affecting the y values that are used during parametric bootstrapping.

OUTPUTS

	b_x_f: Array{Float64}, the bootstrapped x values
	
	b_y_f: Array{Float64}, the bootstrapped y values
	
The bootstrap function can be embedded in a for loop, and will each time produce a different dataset. Performing K times the bootstrapping and fitting each time the model will allow to estimate the error distribution on the peak parameters. This technic has the advantage of making no prior assumption on the probability distribution functions of parameters errors. However, it is  much more time consuming that using the covariance matrix.

The bootperf function
---------------------

This section of the documentation is under work. Here is a quick explanation:

bootperf: reading the bootstrap parameter array for peak fitting

    bootperf(params_boot::Array{Float64}; plotting::ASCIIString = "True", parameter::Int64 = 0, feature::Int64 = 0, histogram_step::Int64 = 100, savefigures::ASCIIString = "False", save_bootrecord::ASCIIString = "Boot_record.pdf", save_histogram::ASCIIString = "Boot_histogram.pdf")

params_boot[i,j] or [i,j,k]: array with i the bootstrrap experiment, j the parameter and k the feature being calculated (e.g., a peak during peak fitting)

plotting: switch to plotting mode ("True") or not ("False"). If true, parameter and feature must be provided, otherwise an error message is returned. 

histogram_step is an integer value to control the histogram X axis division.

savefigures: "True" or "False", explicit, save in the current working directory.

save_bootrecord: Name for the graphic showing the bootstrap mean and std evolutions, with extension.

save_histogram: Name for the graphic showing the histogram for the parameter and feature of interest, with extension.

RETURN: std_record, mean_record, the arrays recording how the standard deviation and mean of the parameters as a function of the bootstrap advance. 

References
----------

For further details, see the following references

Efron, B. 1979. “Bootstrap Methods: Another Look at the Jackknife.” The Annals of Statistics 7 (1): 1–26.

Efron, Bradley. 1981. “Nonparametric Estimates of Standard Error: The Jackknife, the Bootstrap and Other Methods.” Biometrika 68 (3): 589–99. doi:10.1093/biomet/68.3.589.

Efron, B., and Tibshirani, R. 1994. An Introduction to the Bootstrap. CRC press.

	
