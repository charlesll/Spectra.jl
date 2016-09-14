# Spectra.jl

Copyright (c) 2016 Dr. Charles Le Losq

email: charles.lelosq@anu.edu.au

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.53940.svg)](http://dx.doi.org/10.5281/zenodo.53940)

Licence MIT: see LICENCE.md

Spectra.jl is a package aimed at helping spectroscopic (Raman, Infrared, Nuclear Magnetic Resonance, XAS...) data treatment under Julia. Spectra.jl aim is to provide the simplest way to perform actions like baseline fitting and removal, or peak fitting for instance, while respecting the freedom offered by data treatment through using a computer code instead of a Graphic User Interface.

Package is under construction, any help welcome!

# Installation

Installation is done with using the Pkg.add() function of Julia. Please follow the instructions in the documentation (see below). 

The version hosted on Github is bleeding-edge, and will probably NOT WORK well and will change day to day, as I am working on it. I strongly advise not to use it, except if you want to collaborate with me.

# Documentation

Front page for the project is available at http://charlesll.github.io/Spectra.jl/ and full documentation is also available online at http://spectrajl.readthedocs.io/en/latest/ . The latter link provides installation instructions.

# News

Version 0.1.7 is released and contain the following new features:

- The possibility to request different type of normalisation after the correction from the temperature and frequency effects.

- Use of the external calibration mode in RamEau, with added documentation.

Corrections:

- RamEau internal calibration prediction mode equation was wrong and is now corrected


# Quick Usage Examples

A common problem is baseline subtraction and peak fitting when dealing with spectra. After calling the libraries and importing a spectrum like:


	using JuMP
	using PyPlot
	using Ipopt
	using Spectra
	inputsp = readdlm("./examples/data/LS4.txt", '\t') # we import the spectra


We can remove a baseline with the Spectra function "baseline" in two lines:

	roi = [860.0 870.0; 1300.0 1400.0] # the frequencies of the region where we want to fit the baseline
	y_corr, y_bas = baseline(inputsp[:,1],inputsp[:,2],roi,"poly",p=2.0)

A model for fitting gaussian peaks to the spectrum can be easily built with JuMP (https://jump.readthedocs.org/en/latest/):

	mod = Model(solver=IpoptSolver(print_level=0)) # we build a model, initialising the optimiser to Ipopt (https://projects.coin-or.org/Ipopt)
	n = size(x)[1] # number of data
	m = 5 #number of peaks, can be modified!
	@variable(mod,g_amplitudes[i=1:m] >= 0.0) # we declare a first variable containing the amplitudes of peaks
	@variable(mod,g_frequency[i=1:m]) # this one contains the Raman shifts of peaks
	@variable(mod,20.0 <= g_hwhm[i=1:m] <= 40.0) # and this one the hwhm shifts of peaks
	# we set initial values for parameters
	setvalue(g_amplitudes[i=1:m],[1,1,1,1,1])
	setvalue(g_frequency[i=1:m],[950,1050,1090,1140,1190])
	setvalue(g_hwhm[i=1:m],[30,30,30,30,30])
	#We write the model expression that is a sum of Gaussian peaks, and then the objective function that is the least-square deviation function:
	@NLexpression(g_mod[j=1:n],sum{g_amplitudes[i] *exp(-log(2) * ((x[j]-g_frequency[i])/g_hwhm[i])^2), i = 1:m})
	@NLObjective(mod,Min,sum{(g_mod[j] - y[j])^2, j=1:n})
	#And we ask to solve it:
	status = solve(mod)


The "gaussiennes" function allows to get the values fo the peaks after the fit for instance:

	model_peaks, peaks = gaussiennes(amplitudes,frequency,hwhm,x) # we construct the model representation and the individual peaks

Further functions are available, for treating profiles of spectra for instance!

See the examples notebook for futher ideas!

Last modified 14/09/2016
