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

Installation
==================

Two ways of using Spectra.jl: [1] with using a cloud-computing approach and [2] with installing everything on your computer.

[1] JuliaBox (https://www.juliabox.org/) allows you to run Julia in your browser. You still need to add Spectra.jl. To do so, run a notebook, and in the first instance, type

    Pkg.add("Spectra")

Everything shoul install without trouble. Requirements in Spectra.jl are extensive and will provide you all the packages needed by Spectra.jl's functions and examples.

[2] You can download the current version of Julia and follow the installation instruction here: http://julialang.org/downloads/ . Then, run

    Pkg.add("Spectra")

In the Julia shell. Please note that before installing Spectra.jl, the installation of the MatPlotLib library for Python is strongly recommended. Furthermore, some baseline codes call the SciKit learn library, again belonging to the Python ecosystem. Therefore, a good thing will be to install a Python scientific distribution before installing Julia. I recommend Anaconda Python that provides an easy-to-install and nice, fully-featured Python distribution with MatplotLib, SciPy, Numpy and SciKit learn. Follow installation instructions here:

https://www.continuum.io/downloads

## IMPORTANT INFORMATION REGARDING GCVSPLINE ON WINDOWS

Windows users probably need to compile manually the gcvspl.f library if they want to use the GCV spline function. This compilation is automatic on Max OSX and Linux. Please report any problem with that.


Before Starting
==================

Using Julia and Spectra.jl for processing your data is quite similar to Matlab, with the flexibility offered by the open-source and free character of Julia. Reading the docs is strongly recommended. A good start will be to read the docs of Julia itself:

http://docs.julialang.org/en/release-0.4/

Programming can be done locally using your browser and the IJulia notebooks, very similar to the IPython ones. Please follow the instructions there to do so:

https://github.com/JuliaLang/IJulia.jl

For a Matlab-like interface, you can use Atom with Juno. Follow all the instructions at this link to install this interactive developing environment:

http://junolab.org/

For maintaining your packages up-to-date, something critical with the fast evolution of Julia packages, I suggest running each day of Julia use the update command:

    Pkg.update()

For further information about installing packages, please have a read at:

http://docs.julialang.org/en/release-0.4/manual/packages/

Any help developing and maintaining this Spectra.jl package is welcome. You can fork the project on GitHub, modify it and commit your modifications. You can also add requests and everything on Github. Please do not hesitate to do so! The functionalities available in Spectra.jl are not exhaustive, and a little help to add new ones will be more that welcome.

Contents
==================

.. toctree::
   :maxdepth: 2

   PreProcessing
   Integration
   Functions
   Splines
   PeakFitting
   Rameau
   Tutorial
   ToDo
   References

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
