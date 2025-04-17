# Peak fitting

## Introduction

Spectra allows to fit a spectrum with a sum of peaks, of different shapes. For that, we rely on the peak shape functions [`gaussian`](@ref), [`lorentzian`](@ref), [`pseudovoigt`](@ref) and [`pearson7`](@ref). Using those, you can generate a signal given `x`, and peak parameters. To generate a signal composed of multiple contributions from different peaks, use the [`create_peaks`](@ref) function. It is a useful function for instance to create "fake" signals and test our peak fitting function, [`fit_peaks`](@ref). 

[`fit_peaks`](@ref) uses either Optim.jl or a custom quasi-Newton Algorithm to fit the sum of the peaks ``y_{calc}`` to the observed signal ``y``. The regression also takes into account *à priori* errors on model parameters: we assume a Gaussian *prior* on model parameters and this allows constraining the fit. The objective function being minimised is (Tarantola 2005, chapter 3, eq. 3.46):

``\\sum(\\frac{(y - y_{calc})^2}{\\sigma_y^2}) + \\frac{(m_{current} - m_{prior})^2}{\\sigma_m_{prior}^2})``

The least-square regression problem is solved using [BFGS](https://en.wikipedia.org/wiki/Broyden%E2%80%93Fletcher%E2%80%93Goldfarb%E2%80%93Shanno_algorithm) from [Optim.jl](https://julianlsolvers.github.io/Optim.jl/stable/) (`:Optim` backend), or a quasi-Newton (`:qGN` backend) algorithm described in Tarantola ([2005](https://epubs.siam.org/doi/book/10.1137/1.9780898717921), chapter 3). Note that using Optim.jl, you can also set upper and lower boundaries for model parameters. Those are still passed to the function too when using the `:qGN` backend, but in practice they are not used.

[`fit_peaks`](@ref) returns a Named Tuple with various fields, including the areas of the various peaks calculated using [`area_peaks`](@ref).

Errors on model parameters are also available. They are calculated using the inverse of the Hessian matrix at the optimal point, calculated thanks to [ForwardDiff](https://juliadiff.org/ForwardDiff.jl/stable/user/api/#ForwardDiff.hessian). This may not always be the best solution. An alternative way is to calculate the errors with bootstrapping, and for that we have the [`bootsample`] function.

## Peak shapes
```@docs
create_peaks
gaussian
lorentzian
pseudovoigt
pearson7
```

## Fitting function
```@docs
fit_peaks
```

## Utilities
```@docs
area_peaks
bootsample
```

