
# Welcome to Spectra's documentation!

## Introduction

Spectra is a package aimed at helping spectroscopic (Raman, Infrared, Nuclear Magnetic Resonance, XAS...) data treatment written with the [Julia programming language](http://julialang.org/). It's aim is to provide the simplest way to perform actions like baseline fitting and removal or peak fitting for instance, while respecting the freedom offered by data treatment through coding. The key is to provide functions for simplifying the life of the spectroscopist, while still leaving him all the freedom offered by treating data with a performant computer language.

Please consult this documentation to learn using Spectra, do not forget to check the Tips_ section if you have issues, and please report anything you want!

## Starting Notes

Reading the docs is strongly recommended! 

If you are new to Julia, a good start will be to read the [docs of Julia itself](https://docs.julialang.org/en/v1/).

To implement your code, you can modify .jl files in Visual Studio Code for instance and run them from the terminal.

Programming can also be done locally using your browser and the [IJulia Jupyter notebooks](https://github.com/JuliaLang/IJulia.jl), or the reactive [Pluto notebooks](https://plutojl.org/)

## Installation

Two ways of using Spectra: [1] with using a cloud-computing approach and [2] with installing everything on your computer.

[1] [JuliaHub](https://juliahub.com/) allows you to run Julia in your browser. You still need to add Spectra. To do so, run a notebook, and in the first instance, type

```julia
Using Pkg
Pkg.add("Spectra")
```

Everything should install without trouble.

[2] You can install the current version of Julia following the [installation instruction](http://julialang.org/downloads/). Then, type `]` in the Julia REPL shell. You should see `Pkg>` instead of `julia>` on the prompt. If yes, do

```julia-repl
Pkg> add Spectra
```
That's it!

## Getting help

If you wonder about one function, you can use the help mode in Julia's REPL. Just type`?` to enter help mode:
```julia-repl
julia> ?
help?> function_you_want_to_know_about
```

You can also open an [issue on Github](https://github.com/charlesll/Spectra.jl) if you think something is broken, something is missing, an improvement could be made... Do not hesitate! If you do so, please provide details about your setup and if it is about an error, how we can replicate it!

## Contributing

Any help developing and maintaining this Spectra package is welcome. You can fork the project on GitHub, modify it and commit your modifications. You can also add requests and everything on Github. Please do not hesitate to do so! The functionalities available in Spectra are not exhaustive, and a little help to add new ones will be more that welcome.

## Citing Spectra

You can cite Spectra as

	LE LOSQ, C. (2016) Spectra: a Julia package for processing spectroscopic data. Zenodo. 10.5281/zenodo.53940

## Index

The functions that are in Spectra are listed below. See the other part of the documentation for further information.

```@index
```
