**************
Installation
**************

Two ways of using Spectra.jl: [1] with using a cloud-computing approach and [2] with installing everything on your computer.

[1] JuliaBox (https://www.juliabox.org/) allows you to run Julia in your browser. You still need to add Spectra.jl. To do so, run a notebook, and in the first instance, type

    Pkg.add("Spectra")

Everything shoul install without trouble. Requirements in Spectra.jl are extensive and will provide you all the packages needed by Spectra.jl's functions and examples.

[2] You can download the current version of Julia and follow the installation instruction here: http://julialang.org/downloads/ . Then, run

    Pkg.add("Spectra")

In the Julia shell. Please note that before installing Spectra.jl, the installation of the MatPlotLib library for Python is strongly recommended. Furthermore, some baseline codes call the SciKit learn library, again belonging to the Python ecosystem. If not already present in your system, those library should be automatically installed when trying to call for the first time Spectra. However, another good option is to install a Python scientific distribution before installing Julia. I recommend Anaconda Python that provides an easy-to-install and nice, fully-featured Python distribution with MatplotLib, SciPy, Numpy and SciKit learn. Follow installation instructions here:

https://www.continuum.io/downloads

## IMPORTANT INFORMATION REGARDING GCVSPLINE ON WINDOWS

Windows users probably need to compile manually the gcvspl.f library if they want to use the GCV spline function. This compilation is automatic on Max OSX and Linux. Please report any problem with that.

## If you see various errors messages when trying to install Spectra or after a Pkg.update() command, please see the Tips_ section!
