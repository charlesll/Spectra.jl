# SpectraJu

Copyright (c) 2016 Charles Le Losq

email: charles.lelosq@anu.edu.au

Licence MIT: see LICENCE.md

SpectraJu is a package aimed at helping spectroscopic (Raman, Infrared, Nuclear Magnetic Resonance, XAS...) data treatment under Julia.

It's aim is to provide the simplest way to perform actions like baseline fitting and removal or peak fitting for instance, while respecting the freedom offered by data treatment through coding.

It is particularly focused on large datasets because of the high speed of Julia's, e.g. for performing peak fitting along Infrared diffusion profiles.

For peak fitting, the JuMP interface offers a very flexible yet clear way to build models, that can be solve with top-notch solvers such as Ipopt.

Examples are included as notebooks.

Package is under construction, any help welcome!

Last modified 26/02/2016


# Usage Examples

We can start with importing a Raman spectrum (x = Raman shift, y = Intensity) of a glass that we want to peak fit.

We first want to call the libraries we are going to use:
	
	using JuMP
	using PyPlot
	using Ipopt
	using SpectraJu

We import the spectrum and visualize it.

	inputsp = readdlm("./data/LS4.txt", '\t')
	plot(inputsp[:,1],inputsp[:,2],color="black")
	xlabel(L"Raman shift, cm$^{-1}$", fontsize = 14)
	ylabel("Normalized intensity, a. u.", fontsize = 14)
	title("Figure 1: the spectrum of interest",fontsize = 14, fontweight = "bold")

We want to fit the big peak around 1000 cm-1, but first, we need to remove the "background" with the SpectraJu's baseline function.
To do so, we fit a polynomial function at the basis of this peak, and we look at the result:
	roi = [[860.0 870.0]; [1300.0 1400.0]] # the frequencies of the peak basis
	y_corr, y_bas = baseline(inputsp[:,1],inputsp[:,2],roi,"poly",[1.0,1.0,1.0]) # the last vector indicates the coefficients of the baseline k0 + k1 * x + k2 * x^2
	#Creates a plot showing the baseline
	plot(inputsp[:,1],y_corr[:,1],color="blue")
	plot(inputsp[:,1],y_bas[:,1],color="red")
	plot(inputsp[:,1],inputsp[:,2],color="black")
	
We just choose the data between 867 and 1300 cm-1 for peak fitting, and we construct x and y arrays accordingly (TODO: write a function for that):
	index_interest = find(867.0 .< inputsp[:,1] .< 1300.0)
	interestspectra = y_corr[index_interest,1]
	ese0 = sqrt(abs(interestspectra[:,1]))/abs(interestspectra[:,1]) # the relative errors after baseline subtraction
	interestspectra[:,1] = interestspectra[:,1]/trapz(inputsp[index_interest,1],interestspectra[:,1])*100 # normalise spectra to maximum intensity, easier to handle 
	# First we simplify things by calling x, y and the frequency and intensity of spectra for later use
	sigma = abs(ese0.*interestspectra[:,1]) #calculate good ese
	x = inputsp[index_interest,1]
	y = interestspectra[:,1]
	println("Done")

Then we construct the model for fitting gaussian peaks to the spectrum using JuMP (https://jump.readthedocs.org/en/latest/)
	mod = Model(solver=IpoptSolver(print_level=0)) # we build a model, initialising the optimiser to Ipopt (https://projects.coin-or.org/Ipopt)
	n = size(x)[1] # number of data
	m = 5 #number of peaks, can be modified!

	@defVar(mod,g_amplitudes[i=1:m] >= 0.0) # we declare a first variable containing the amplitudes of peaks
	@defVar(mod,g_frequency[i=1:m]) # this one contains the Raman shifts of peaks
	@defVar(mod,20.0 <= g_hwhm[i=1:m] <= 40.0) # and this one the hwhm shifts of peaks

You notice in the last part, we also added little constrains for the amplitudes and the hwhm. Those are soft constrains, but necessary to have good results for this kind of "broad band" fit. Next, we initiate the values to our guess. The more important thing here is the frequencies, then the hwhm, and finally the amplitudes (not critical).

	setValue(g_amplitudes[i=1:m],[1,1,1,1,1])
	setValue(g_frequency[i=1:m],[950,1050,1090,1140,1190])
	setValue(g_hwhm[i=1:m],[30,30,30,30,30])

We write the model expression that is a sum of Gaussian peaks, and then the objective function that is the least-square deviation function:
	@defNLExpr(g_mod[j=1:n],sum{g_amplitudes[i] *exp(-log(2) * ((x[j]-g_frequency[i])/g_hwhm[i])^2), i = 1:m})
	@setNLObjective(mod,Min,sum{(g_mod[j] - y[j])^2, j=1:n}) 

And we ask to solve it:
	status = solve(mod)
	
We get the values of adjusted peak parameters and plot the result, with using the SpectraJu gaussiennes function to get the peaks:
	# parameter extractions
	amplitudes = getValue(g_amplitudes)
	frequency = getValue(g_frequency)
	hwhm = getValue(g_hwhm)

	model_peaks, peaks = gaussiennes(amplitudes,frequency,hwhm,x) # we construct the model representation and the individual peaks

	#we plot the results
	plot(x,y,color="blue")
	plot(x,model_peaks,color="blue")
	plot(x,peaks[:,:],color="red")

And voila! To be implemented in a very close future: error estimation with bootstrapping the data!
