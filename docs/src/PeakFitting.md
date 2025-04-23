# Peak fitting

## Introduction

Spectra allows to fit a spectrum with a sum of peaks, of different shapes. For that, we rely on the peak shape functions [`gaussian`](@ref), [`lorentzian`](@ref), [`pseudovoigt`](@ref) and [`pearson7`](@ref). Using those, you can generate a signal given `x`, and peak parameters. To generate a signal composed of multiple contributions from different peaks, use the [`create_peaks`](@ref) function. It is a useful function for instance to create "fake" signals and test our peak fitting function, [`fit_peaks`](@ref). 

[`fit_peaks`](@ref) uses either Optim.jl or a custom quasi-Newton algorithm to fit the sum of the peaks `` y_{calc} `` to an observed signal `` y `` affected by errors `` \sigma_y ``. The regression also takes into account *à priori* errors `` \sigma_{m_{prior}} `` on prior model parameters `` m_{prior} ``: we assume a Gaussian *prior* on model parameters with a mean ``m_{prior}`` and a covariance matrix ``C_M`` which diagonal contains ``\sigma_{m_{prior}}^2``.

Given the forward calculation of `` y_{calc} `` as

`` y_{calc} = g(m) ``

with ``m`` the model parameters, the misfit function ``S`` is ([Tarantola 2005](https://epubs.siam.org/doi/book/10.1137/1.9780898717921), chapter 3, eq. 3.46):

`` S(m) = \frac{1}{2}[(y - y_{calc})^{t} C_D^{-1}(y - y_{calc}) + (m - m_{prior})^{t}C_M^{-1}(m-m_{prior})] ``

where ``C_D^{-1}`` is the inverse of the data covariance matrix, which diagonal contains ``\sigma_y^2``. The misfit function is related to the posterior probability density in the model space following ``K \exp (-S(m))``, with ``K`` a constant. Here, we are dealing with a non-linear problem so this posterior probability density is not Gaussian. However, we assume that it can be linearized near ``m_{prior}``. Therefore, starting close to or at ``m_{prior}``, we can use the following quasi-Newton algorithm to find a suitable solution (Tarantola 2005, eq. 3.51):

`` m_{n+1} = m_{n} - \mu_n(G_n^tC_D^{-1}G_n + C_M^{-1})^{-1}(G_n^tC_D^{-1}(y_{calc, n} - y) + C_M^{-1}(m_n - m_{prior})) ``

where `` y_{calc, n} `` is the model output at iteration `` n ``, `` (G_n)^i_\alpha = (\frac{g^i}{m^\alpha})_{m_n} `` is the matrix of partial derivatives, and ``\mu_n`` a step size typically lower or equal to 1. In `fit_peaks`, the two parameters `maxiter` and `relax` control the maximum number of iterations `n` and the step size `` \mu_n ``, with `` \mu_n = \frac{1}{\\text{relax}} ``.

The other available method is the Interior Point Newton algorithm from [Optim.jl](https://julianlsolvers.github.io/Optim.jl/stable/algo/ipnewton/) (`:Optim` backend). This method implements an interior-point primal-dual Newton algorithm for solving the nonlinear, constrained optimization problem. In other terms, it allows the use of constraints, such as parameter boundaries. This is the main difference with the quasi-Newton method. The misfit function that is minimized is the ``S(m)`` function provided above.

Convergence between the two methods usually is similar, with the quasi-Newton method being slightly faster. However, as there is no boundaries on parameters, one may have problems for instance if the intensity of one peak is close to 0. The quasi-Newton method offers no boundaries for the parameter values, contrary to the Interior Point Newton method. 

## Fitting procedure

A short example probably is better than many worlds. First, you define a vector containing named vectors of peak types, ``m_{prior}``, ``\sigma_{m_{prior}}``, and lower and upper boundaries. For instance, we declare a vector of peak informations like this:

```julia 
peaks_info = [
        ## (type, initial_params, uncertainties, lower_bounds, upper_bounds)
        (:gaussian, # peak type   
        [0.20, 950, 27], # m_prior (intensity, position, hwhm)
        [0.05, 5.0,3.0], # sigma_m_prior
        [0.0, 0.0, 0.0], # lower boundaries
        [Inf, Inf, 60.0]), # upper boundaries 

        (:lorentzian, # peak type        
        [0.44, 1044, 40], # m_prior (intensity, position, hwhm)
        [0.05, 5.0, 3.0], # sigma_m_prior
        [0.0, 0.0, 0.0], # lower boundaries
        [Inf, Inf, 60.0]), # upper boundaries

        (:pseudovoigt, # peak type  
        [0.80, 1086, 30., 0.8], # m_prior (intensity, position, hwhm, lorentzian_ratio)
        [0.05, 5.0, 3.0, 0.02], # sigma_m_prior
        [0.0, 0.0, 0.0, 0.0], # lower boundaries
        [Inf, Inf, 60.0, 1.0]), # upper boundaries

        (:pearson7, # peak type  
        [0.50, 1086, 30., 0.8], # m_prior (intensity, position, hwhm, shape exponent)
        [0.01, 5.0, 3.0, 0.02], # sigma_m_prior
        [0.0, 0.0, 0.0, 0.0], # lower boundaries
        [Inf, Inf, 60.0, 1.0]), # upper boundaries
    ]
```
!!! note

    Upper and lower boundaries for model parameters are always required when constructing this peaks information vector. However, remember that they are not used in the quasi-Newton method.

We then pass the data and this vector of peak informations to [`prepare_context`](@ref) that stores everything in a Julia object. Given `x`, `y` and `sigma_y` the data vectors, we write

```julia
ctx = prepare_context(x, y, peaks_info, sigma_y)
```

!!! note

    `sigma_y` is an optional argument that will be set to an array of 1 if you do not pass a vector. It is advised to pass proper errors for a proper scaling of the misfit function.

From there, a good thing is to check that your prior model is not too remote from a good fit. The algorithms we use are local optimization methods and do not aim at finding a global minimum, but only local solutions. If you set ``m_{prior}`` to values far from the solution, they will fail. To check your starting parameters, we can fit the prior model and the data using:

```julia
p = plot_fit(ctx, title="Prior model")
display(p)
```

Modify your starting parameters until the model starts to make sense, and then perform a fit calling [`fit_peaks`](@ref) with your favorite backend. For instance, if you want to use IPNewton, call:

```julia
result = fit_peaks(ctx, backend=:Optim)
```

or if you want to use the quasi-Newton method described above, call:


```julia
result = fit_peaks(ctx, backend=:qGN, maxiter=100, relax=5)
```

`result` is an object containing:
- `context::FitContext`: Fit context
- `peak_results::Vector`: Peak results
- `params::Vector{Float64}`: Peak parameters with uncertainties
- `covariance::Matrix{Float64}`: Covariance matrix
- `errors::Vector{Float64}`: 1-sigma standard errors on parameters
- `y_calc::Vector{Float64}`: Model predictions
- `residuals::Vector{Float64}`: Residuals

In `result.peal_results`, numbers are now Measurements as we use the [Measurements.jl](https://github.com/JuliaPhysics/Measurements.jl) library to automatically propagate fitting errors for peak area calculations. 

You can print the parameters and peak areas using [`print_params`](@ref) and plot the fit using [`plot_fit`](@ref): 

```julia
print_params(result.peak_results)

plot_fit(ctx, result.peak_results)
```

## Errors on peak parameters

###  Errors provided by `fit_peaks`

The errors provided by `fit_peaks` come from the evaluation of the Hessian matrix at the optimal point. In the quasi-Newton algorithm, we directly calculate the posterior model covariance matrix as (Tarantola, 2005, eq. 3.53):

`` C_{M, post} = C_{M} - C_{M}G^{t}(G C_M G{t} +C_D)^{-1}G C_M ``

and retrieve the standard errors on model parameters from the squared root of the diagonal of ``C_{M, post}``. In the IPNewton case, we use ForwardDiff.hessian() to calculate the Hessian matrix at the optimal point to optain ``C_{M, post}``.

!!! note

    Do not neglect off-diagonal terms in ``C_{M, post}``. Peak parameters often are strongly correlated, and neglecting off-diagonal terms may lead to report wrong errors on quantities that depends on several peak parameters, such as peak areas. This is why we automatically propagate errors using ``C_{M, post}`` and Measurements.jl when calculating peak areas in `fit_peaks` for instance.

### Checking errors with bootstrapping

``C_{M, post}`` may not necessarily contain valid parameter uncertainties, particularly if the problem is strongly non-linear.

To check for parameter errors, one option is to use [bootstrapping](https://en.wikipedia.org/wiki/Bootstrapping_(statistics)).
[`bootstrap`](@ref) allows us to bootstrap a spectrum and refit the model on the spectrum subsamples. You will obtain samples of models all adjusted on slightly different data. Using the model samples, we can check that the errors calculated from the Hessian are valid, or if they are not, we can use the new errors calculated from the bootstrap samples.

The interface is easy, similar to that of `fit_peaks` but the fit context is defined internally so 
you don't even have to worry about that. For instance, using the quasi-Newton method, we can do: 

```julia
boot_params, boot_results = bootstrap(x, y, sigma_y, peaks_info, nb_boot = 1000, backend=:qGN, relax=5., maxiter=20)
```
!!! tip

    The quasi-Newton algorithm is the quickest so you may prefer this one for bootstrapping. If so, you can try also setting maxiter to a value as small as possible without affecting fit convergence. Another tip is to provide a new `peaks_info` vector with ``m_{prior}`` set to the values of a good fit that you previously did.

The `bootstrap` function returns:
  - a matrix of size (nb_params, nb_boot) with the fitted parameters (here `boot_params`);
  - a peak_results objects with values tied to their errors thanks to Measurements.jl (here called `boot_results`).

We can now print the bootstrapped results and compare the errors with those previously calculated from the Hessian:

```julia
print_params(boot_results)
```
OK, actually for this example, we see that the errors from the boostrap analysis
are close to those calculated from the Hessian matrix. Everything thus seems OK.

Of course, a final quick visual inspection is always nice. This can be done by passing the 
median of the matrix of bootstrapped parameters to the plot_fit function:

plot_fit(ctx, boot_results)

## Final remarks

Done, please do check the examples in the Tutorials section for further peak fitting examples. Below you will find the full API of the various functions, including peak shapes, area calculations, fitting algorithms...

## Peak fitting functions

```@docs
prepare_context
fit_peaks
fit_Optim
fit_qNewton
plot_fit
print_params
get_peak_results
bootstrap
FitContext
FitResult
```

## Peak shapes
```@docs
create_peaks
gaussian
lorentzian
pseudovoigt
pearson7
```

## Utilities
```@docs
area_peaks
bootsample
```

