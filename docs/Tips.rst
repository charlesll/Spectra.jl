.. _Tips:
***********************
Tips
***********************

In this section are listed various tips for the use of Julia and Spectra:

------------
Installation
------------

1) You need the gfortran, gcc and g++ compilers. ifort also works. Check that you have them on your system. Even Ubuntu does not necessary come with those compilers out of the box. If you don't know anything about installing them, ask Google: "Installing gcc/gfortran/g++ on my mac/linux/windows"

2) Windows users probably need to manually compile the gcvspl.f library if they want to use the GCV spline function. This compilation is automatic on Max OSX and Linux. Please report any problem with that.

3) If you see errors messages linked to PyCall, you may have a problem with your environment variable. To solve it, tyope the following commands in the Julia prompt:

	ENV["PYTHON"]=""
	Pkg.build("PyCall")

At this point it should work. If yes, you now can enter:

	Pkg.add("Spectra")

---------------
Maintenance
---------------

1) The Julia package ecosystem is constantly evolving, with daily changes. Because of that, it is strongly recommanded to run in the starting Julia prompt a

	Pkg.update()

command every day, at the beginning of your session.

2) Time to time, after running the Pkg.update() command for instance and trying to directly work with the same Julia session, you may get Warning/Error messages during the packages pre-compilation indicating a problem with Compat. To solve that, just quit the current session (exit the notebooks AND close the terminals), and open a new Julia terminal. Most of the time, this solves the problem.

---------------
Running Spectra
---------------

1) Spectra is changing every week, if not every day in some case. Do not forget to Pkg.update() quite often, and check the website.

2) Always be careful to enter float and integer numbers as required by the functions! They will return an error message if you do not do that.

3) For the spline, do not hesitate to test a broad range in term of order of magnitudes for the smoothing parameter.

4) SVMregression and KRregression will take more time as several models are tried over a broad range of hyperparameters. Therefore, it is normal that those technics require more time, up to ten to twenty minutes for treating 50 to 100 spectra.

------------------
Potential problems
------------------

1) Using Julia on a Fedora Linux installed in a VirtualBox virtual machine, I encountered the issue of memory mapping not working when trying to read with `readdlm`/`readcsv` some files that where in a VirtualBox shared folder:

```
LoadError: SystemError: memory mapping failed: Invalid argument
```

This issue is solved by setting the optinal argument use_mmap = false in the readcsv/readdlm call. There is an option `mmap_switch` (false/true) in `rameau` that allows also to set use_mmap, in case you encouter this problem when calling `rameau` in a virtual environment.
