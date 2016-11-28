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

For Windows users, Spectra.jl will issue a WARNING message saying that GCVSPL.F is not compiled automatically upon installation, and will point to this page. You will need to compile GCVSPL.F by yourself for now. If you want to avoid this step, I recommand using JuliaBox.org where everything can run smoothly, or using Julia inside a free virtualbox Linux installation (https://www.virtualbox.org/). This makes things pretty easy. If you want to run Julia directly on your Windows system, you can try the following steps to compile GCVSPL.F with cygwin:

	1) create bin32 and bin64 folders in the /deps forlder;
	
	2) compile GCVSPL.F as a shared libgcvspl.dll library in ./bin32 or ./bin64. Using cygwin, this can be done as:
	
	i686-w64-mingw32-gfortran -o bin32/libgcvspl.dll -O3 -shared -static-libgfortran -static-libgcc src/gcvspline/*.f
	
	x86_64-w64-mingw32-gfortran -o bin64/libgcvspl.dll -O3 -shared -static-libgfortran -static-libgcc src/gcvspline/*.f
	
	3) if this is not working, you may want to also change the winpath in Spectra.jl,  see /Spectra/src/Spectra.jl line 38.
	
I never tested those steps because I do not have a Windows system available, so I am not sure if they fully work. You might have to tweak things a little bit. This will be corrected soon. If anybody would like to help me with that, please submit a pull request of a working Windows installation procedure.

## If you see various errors messages when trying to install Spectra or after a Pkg.update() command, please see the Tips_ section!
