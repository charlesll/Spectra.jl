# Installation

## General Instructions

Two ways of using Spectra: [1] with using a cloud-computing approach and [2] with installing everything on your computer.

[1] JuliaBox (https://www.juliabox.org/) allows you to run Julia in your browser. You still need to add Spectra. To do so, run a notebook, and in the first instance, type

	```
	Using Pkg
	Pkg.add("Spectra")
	```

Everything should install without trouble. Requirements in Spectra are extensive and will provide you all the packages needed by Spectra's functions and examples.

[2] You can download the current version of Julia and follow the installation instruction here: http://julialang.org/downloads/ . Then, type `]` in the Julia REPL shell. You should see `Pkg>` instead of `julia>` on the prompt. If yes, do

    ```julia-repl
    Pkg> add Spectra
    ```

### Notes on gcvspline for v0.4.1 and higher

Installing the gcvspline Python library is optional but gives access to the GCV algorithms in the smooth() and baseline() functions. Install as

    ```julia-repl
     julia> using PyCall
     julia> pyimport_conda("pip", "pip")
     julia> run(`$(PyCall.python) -m pip install $("gcvspline")`)
     ```

#### Prior to v0.4.1:

Please update to version v0.4.1 :)

## Error messages?

If you see various errors messages when trying to install Spectra or after a Pkg.update() command, please see the Tips section!
