
# Welcome to Spectra's documentation!

## Introduction

Spectra.jl is a Julia package designed to simplify the treatment of spectroscopic data (Raman, Infrared, Nuclear Magnetic Resonance, XAS, etc.). It provides straightforward tools for common tasks such as smoothing, baseline fitting, and peak fitting, while preserving the flexibility and power of programmatic data analysis in Julia. Spectra.jl is particularly suited for large datasets and integrates well with other Julia packages, such as JuMP for model building.

Please consult this documentation to learn how to use Spectra.jl. Check the [Tips](@ref) section if you encounter issues, and feel free to report problems or suggest improvements!

## Starting Notes

- Reading the docs is strongly recommended! 

- If you are new to Julia, start by reading the [official documentation of Julia](https://docs.julialang.org/en/v1/).

- You can develop your code in text editors such as Visual Studio Code, or use interactive environments such as [IJulia Jupyter notebooks](https://github.com/JuliaLang/IJulia.jl) or the reactive [Pluto notebooks](https://plutojl.org/).

!!! important

	Spectra v2.0.0 API has evolved a lot since v1.0.0. There are breaking changes that may require to update your codes.

## Installation

**Cloud (JuliaHub):**
1. Launch a notebook on [JuliaHub](https://juliahub.com/).
2. In your first cell, run:
    ```
    using Pkg
    Pkg.add("Spectra")
    ```

**Local:**
1. Download Julia from the [official website](http://julialang.org/downloads/).
2. Open the Julia REPL, press `]` to enter package mode (`Pkg>`), then run:
    ```
    Pkg> add Spectra
    ```

## Getting help

- The search bar on the top left of this website is useful to search the documentation, or also you can use the function index at the end of this page for instance.

- You can also use the help mode in Julia's REPL by typing `?` in the REPL, then the function name:

```julia-repl
julia> ?
help?> function_you_want_to_know_about
```

- If you encounter issues or have suggestions, open an [issue on GitHub](https://github.com/charlesll/Spectra.jl) with details about your setup and steps to reproduce any errors.

## Contributing

Contributions are welcome! You can fork the project, submit pull requests, or open issues for bugs and feature requests. The feature set is not exhaustive, and community help is appreciated.

## Citing Spectra

If you use Spectra.jl in your work, please cite:

	LE LOSQ, C. (2016) Spectra: a Julia package for processing spectroscopic data. Zenodo. 10.5281/zenodo.53940

## Index

The functions available in Spectra.jl are listed below. See the rest of the documentation for detailed usage.

```@index
```
