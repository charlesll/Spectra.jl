# Spectra.jl

Copyright (c) 2016 Charles Le Losq

email: charles.lelosq@anu.edu.au

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.53940.svg)](http://dx.doi.org/10.5281/zenodo.53940)

Licence MIT: see LICENCE.md

Spectra.jl is a package aimed at helping spectroscopic (Raman, Infrared, Nuclear Magnetic Resonance, XAS...) data treatment under Julia.

See the github wiki for further information.

Spectra.jl aim is to provide the simplest way to perform actions like baseline fitting and removal, or peak fitting for instance, while respecting the freedom offered by data treatment through using a computer code instead of a Graphic User Interface.

It is particularly focused on large datasets because of the high speed of Julia's, e.g. for performing peak fitting along Infrared diffusion profiles.

For peak fitting, the JuMP interface offers a very flexible yet clear way to build models, that can be solve with top-notch solvers such as Ipopt.

Examples are included as notebooks.

Package is under construction, any help welcome!

Last modified 19/03/2016


# Usage Examples

A common problem is baseline subtraction and peak fitting when dealing with spectra. After calling the libraries and importing a spectrum like:


	using JuMP
	using PyPlot
	using Ipopt
	using Spectra
	inputsp = readdlm("./examples/data/LS4.txt", '\t') # we import the spectra


We can remove a baseline with the Spectra function "baseline" in two lines:

	roi = [[860.0 870.0]; [1300.0 1400.0]] # the frequencies of the region where we want to fit the baseline
	y_corr, y_bas = baseline(inputsp[:,1],inputsp[:,2],roi,"poly",[1.0,1.0,1.0]) # the last vector indicates the coefficients of the baseline k0 + k1 * x + k2 * x^2

A model for fitting gaussian peaks to the spectrum can be easily built with JuMP (https://jump.readthedocs.org/en/latest/):

	mod = Model(solver=IpoptSolver(print_level=0)) # we build a model, initialising the optimiser to Ipopt (https://projects.coin-or.org/Ipopt)
	n = size(x)[1] # number of data
	m = 5 #number of peaks, can be modified!
	@defVar(mod,g_amplitudes[i=1:m] >= 0.0) # we declare a first variable containing the amplitudes of peaks
	@defVar(mod,g_frequency[i=1:m]) # this one contains the Raman shifts of peaks
	@defVar(mod,20.0 <= g_hwhm[i=1:m] <= 40.0) # and this one the hwhm shifts of peaks
	# we set initial values for parameters
	setValue(g_amplitudes[i=1:m],[1,1,1,1,1])
	setValue(g_frequency[i=1:m],[950,1050,1090,1140,1190])
	setValue(g_hwhm[i=1:m],[30,30,30,30,30])
	#We write the model expression that is a sum of Gaussian peaks, and then the objective function that is the least-square deviation function:
	@defNLExpr(g_mod[j=1:n],sum{g_amplitudes[i] *exp(-log(2) * ((x[j]-g_frequency[i])/g_hwhm[i])^2), i = 1:m})
	@setNLObjective(mod,Min,sum{(g_mod[j] - y[j])^2, j=1:n})
	#And we ask to solve it:
	status = solve(mod)


The "gaussiennes" function allows to get the values fo the peaks after the fit for instance:

	model_peaks, peaks = gaussiennes(amplitudes,frequency,hwhm,x) # we construct the model representation and the individual peaks

Further functions are available, for treating profiles of spectra for instance!

See the examples notebook for futher ideas!
