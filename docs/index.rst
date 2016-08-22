.. Spectra.jl documentation master file, created by
   sphinx-quickstart on Wed Jun  1 10:49:17 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Spectra.jl's documentation!
======================================

Introduction
==================
Spectra.jl is a package aimed at helping spectroscopic (Raman, Infrared, Nuclear Magnetic Resonance, XAS...) data treatment written with the Julia programming language (http://julialang.org/). It's aim is to provide the simplest way to perform actions like baseline fitting and removal or peak fitting for instance, while respecting the freedom offered by data treatment through coding. Therefore, Spectra.jl is aimed to be used explicitly with other packages like JuMP for building models (see http://www.juliaopt.org/ for further details on optimisation under Julia). The key is to provide functions for simplifying the life of the spectroscopist, while still leaving him all the freedom offered by treating data with a performant computer language.

Spectra.jl is particularly focused on large datasets because of the high speed of Julia's, e.g. for performing peak fitting along Infrared diffusion profiles. For peak fitting for instance, the JuMP interface offers a very flexible yet clear way to build models, that can be solve with solvers such as Ipopt or NLopt.

Please consult this documentation to learn using Spectra, do not forget to check the Tips_ section if you have issues, and please report anything you want!

Contents
==================

.. toctree::
   :maxdepth: 2

   Beforestart
   Installation
   Tips
   Tutorial
   PreProcessing
   Integration
   Functions
   Splines
   PeakFitting
   Rameau
   ToDo
   References

Citing Spectra
==============

You can cite Spectra as 

LE LOSQ, C. (2016) Spectra.jl: a Julia package for processing spectroscopic data. Zenodo. 10.5281/zenodo.53940


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
