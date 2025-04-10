
# Welcome to Spectra's documentation!

## Introduction

Spectra is a package aimed at helping spectroscopic (Raman, Infrared, Nuclear Magnetic Resonance, XAS...) data treatment written with the [Julia programming language](http://julialang.org/). It's aim is to provide the simplest way to perform actions like baseline fitting and removal or peak fitting for instance, while respecting the freedom offered by data treatment through coding. Therefore, Spectra is aimed to be used explicitly with other packages like [JuMP](http://www.juliaopt.org/) for building models. The key is to provide functions for simplifying the life of the spectroscopist, while still leaving him all the freedom offered by treating data with a performant computer language.

Spectra is particularly focused on large datasets because of the high speed of Julia's, e.g. for performing peak fitting along Infrared diffusion profiles. For peak fitting for instance, the JuMP interface offers a very flexible yet clear way to build models, that can be solve with solvers such as Ipopt or NLopt.

Please consult this documentation to learn using Spectra, do not forget to check the Tips_ section if you have issues, and please report anything you want!

## Starting Notes

Using Julia and Spectra for processing your data is quite similar to Matlab, with the flexibility offered by the open-source and free character of Julia. Reading the docs is strongly recommended. A good start will be to read the [docs of Julia itself](http://docs.julialang.org/en/release-0.5/).

Programming can be done locally using your browser and the [IJulia notebooks](https://github.com/JuliaLang/IJulia.jl), very similar to the IPython ones. For a Matlab-like interface, you can use [Atom with Juno](http://junolab.org/).

For maintaining your packages up-to-date, something critical with the fast evolution of Julia packages, I suggest running each day of Julia use the update command:

	Pkg.update()

## Installation

Two ways of using Spectra: [1] with using a cloud-computing approach and [2] with installing everything on your computer.

[1] [JuliaHub](https://juliahub.com/) allows you to run Julia in your browser. You still need to add Spectra. To do so, run a notebook, and in the first instance, type

```julia
Using Pkg
Pkg.add("Spectra")
```

Everything should install without trouble. Requirements in Spectra are extensive and will provide you all the packages needed by Spectra's functions and examples.

[2] You can download the current version of Julia and follow the [installation instruction](http://julialang.org/downloads/). Then, type `]` in the Julia REPL shell. You should see `Pkg>` instead of `julia>` on the prompt. If yes, do

```julia-repl
Pkg> add Spectra
```

## Examples

Examples are available in the [examples folder of Spectra](https://github.com/charlesll/Spectra/tree/master/examples).

## Contributing

Any help developing and maintaining this Spectra package is welcome. You can fork the project on GitHub, modify it and commit your modifications. You can also add requests and everything on Github. Please do not hesitate to do so! The functionalities available in Spectra are not exhaustive, and a little help to add new ones will be more that welcome.

## Citing Spectra

You can cite Spectra as

	LE LOSQ, C. (2016) Spectra: a Julia package for processing spectroscopic data. Zenodo. 10.5281/zenodo.53940

## Index

The functions that are in Spectra are listed below. See the other part of the documentation for further information.

```@index
```
