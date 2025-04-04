# Installation

## General Instructions

Two ways of using Spectra: [1] with using a cloud-computing approach and [2] with installing everything on your computer.

[1] [JuliaHub](https://juliahub.com/) allows you to run Julia in your browser. You still need to add Spectra. To do so, run a notebook, and in the first instance, type

	```
	Using Pkg
	Pkg.add("Spectra")
	```

Everything should install without trouble. Requirements in Spectra are extensive and will provide you all the packages needed by Spectra's functions and examples.

[2] You can download the current version of Julia and follow the [installation instruction](http://julialang.org/downloads/). Then, type `]` in the Julia REPL shell. You should see `Pkg>` instead of `julia>` on the prompt. If yes, do

    ```julia-repl
    Pkg> add Spectra
    ```

