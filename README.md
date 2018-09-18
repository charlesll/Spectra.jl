# Spectra.jl

Copyright (c) 2016-2017 Dr. Charles Le Losq

email: charles.lelosq@anu.edu.au

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.53940.svg)](http://dx.doi.org/10.5281/zenodo.53940)

Licence MIT: see LICENCE.md

Spectra.jl is a package aimed at helping spectroscopic (Raman, Infrared, Nuclear Magnetic Resonance, XAS...) data treatment under Julia. Spectra.jl aim is to provide the simplest way to perform actions like baseline fitting and removal, or peak fitting for instance, while respecting the freedom offered by data treatment through using a computer code instead of a Graphic User Interface.

Package is under construction, any help welcome!

# Installation

Installation is done with using the Pkg.add() function of Julia. Please follow the instructions in the documentation (see below). 

The version hosted on Github is bleeding-edge, and will probably NOT WORK well and will change day to day, as I am working on it. I strongly advise not to use it, except if you want to collaborate with me.

# Documentation

For old versions of Spectra <= 0.3.1, please see the documentation http://spectrajl.readthedocs.io/en/v0.3.1/ and earlier related versions.

Full documentation for the stable and latest versions is now available at

[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://charlesll.github.io/Spectra.jl/latest) [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://charlesll.github.io/Spectra.jl/stable)

# News

See the NEWS.md file for following the evolution of Spectra.jl.

# Installation

In Julia 0.6 or earlier, use `Pkg.add("Spectra")`.

In Julia 0.7 and later, use the `pkg` environment (key `]`), then directly run `add git@github.com:charlesll/Spectra.jl.git`. For those versions of Julia, you need to install through this way as not all dependencies have solved version problems and there is many warnings in Spectra, such that I can't tag a proper 0.4.0 version of Spectra yet. You also will have to manually install NMF using `add NMF#aa/1.0`.

# Quick Usage Examples

A common problem is baseline subtraction and peak fitting when dealing with spectra. After calling the libraries and importing a spectrum like:


	using JuMP, PyPlot, Spectra

	data = readdlm("./examples/data/LS4.txt", '\t') # we import the spectra
	x = data[:,1]
	y = data[:,2]


We can remove a baseline with the Spectra function "baseline" in two lines:

	roi = [860.0 870.0; 1300.0 1400.0] # the frequencies of the region where we want to fit the baseline
	y_corr, y_bas = baseline(x,y,roi,"poly",p=2.0)

Further functions are available, for treating profiles of spectra for instance!

See the examples notebook for futher ideas!

Last modified May 2017
