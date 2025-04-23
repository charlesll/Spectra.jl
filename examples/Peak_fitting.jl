# # Peak fitting
# 
# Charles Le Losq, May 2017; Updated January 2024, April 2025
# 
# In the new version, we showcase the use of the `fit_peaks` function and then show also how you can solve the problem with Turing.jl and JuMP too. Running the notebook entirely thus requires you to have those installed, but if you want to stop at the `fit_peaks` function, you can skip the cells related to Turing.jl and JuMP!
# 
# ## Generating data
# 
# We generate three peaks using Spectra's functions:
# 
using Plots, Spectra, StatsBase

## The X axis
x = collect(0:0.2:100)
## The "perfect" Y signal
y_perfect = (
    gaussian(x, [10.0, 35.0, 10.0]) +
    lorentzian(x, [15.0, 55.0, 3.0]) +
    pearson7(x, [20.0, 45.0, 2.0, 0.4])
)
## Of course, in real world we have noise, here it is Gaussian
noise = randn(length(x))*0.3

## This is what we observe and want to fit
y_obs = y_perfect + noise

## Let's visualize it!
p1 = plot(x, [y_perfect, y_obs]; labels=["The true signal" "Observations"])
savefig("fit_1.svg"); nothing #hide
# ![](fit_1.svg)

# ## Declaring peak informations
#
# In the real world, we see only the above figure. 
# We guess peak positions and shapes based on the visual inspection. 
# So we start with prior estimates on peak parameters based on this visual inspection. 
# Data from publications can bring further insights. This helps e.g. constraining the position of 
# some peaks, and `Spectra.fit_peaks` leverages such information via an *à priori* vector of 
# peak parameter uncertainties (1-sigma).
# 
# Another option is available : you can place lower and upper bounds on parameter values. 
# This is very useful if your parameters cannot take some values: e.g. hwhm should be positive 
# and lorentzian_fraction should be comprised between 0 and 1.
# 
# All this is indicated as follow:

peaks_info = [
    ## (peak_type, initial_params, params_uncertainties, lower_bounds, upper_bounds)
    (   
        :gaussian, ## peak_type
        [10.5, 30.0, 11.0], ## initial_params > prior model
        [5.0, 5.0, 3.0], ## params_uncertainties  > prior model uncertainties
        [0.0, 0.0, 0.0], ## lower_bounds
        [Inf, Inf, 50.0]), ## upper_bounds
    (
        :lorentzian, ## peak_type
        [17.5, 54.0, 3.1], ## initial_params > prior model
        [5.0, 3.0, 1.0], ## params_uncertainties > prior model uncertainties
        [0.0, 0.0, 0.0], ## lower_bounds
        [Inf, Inf, Inf], ## upper_bounds
    ),
    (
        :pearson7, ## peak_type
        [21.5, 44.0, 3.0, 0.4], ## initial_params > prior model
        [3.0, 2.0, 5.0, 0.02], ## params_uncertainties > prior model uncertainties
        [0.0, 0.0, 0.0, 0.0], ## lower_bounds
        [100.0, 100.0, 50.0, Inf], ## upper_bounds
    ),
]

# ## Error on the data?
#
# Usually we also do not have a precise idea of the errors that affect y. A good approximation 
# can be to simply smooth the signal, and calculate the errors using the smoothed signal and the observed one.
# Let's also do that here and check that it yields a fairly good approximation of the true error (should be around 0.3)

y_smo = smooth(x, y_obs, method="gcvspline");
estimated_mean_error = sqrt(mean((y_obs .- y_smo).^2))
println("The estimated mean standard error on y_obs is $(round(estimated_mean_error,digits=2))")

# OK, the result seems to be not too bad. We will place ourselves in a "real world" situation and 
# use those errors. For convenience we create an vector of errors

estimated_error = estimated_mean_error * ones(size(x));

# ## Fitting
#
# We use now the `fit_peaks` function to perform the fit.
# 
# This function returns a Named Tuple with various fields (see documentation).
# 
# First, we need to create the Julia object that contains all the context of the fit, 
# including the data, priors, peak informations...

ctx = prepare_context(x, y_obs, peaks_info, estimated_error)

# We launch the fit using the quasi-Newton algorithm implemented in Spectra. 
# This algorithm assumes that everything is Gaussian: the distributions
# of parameter and data errors are assumed to be Gaussian, and we also
# assume that the problem is nearly linear near the optimal point. In practice, 
# this algorithm works well and has the advantage of leveraging the prior information provided
# in `peaks_info`.
# The downside is sometimes some numerical instability: when `relax`is too low, DomainError is raised because 
# complex numbers appear in the Hessian matrix and mess with the sqrt() function. The fit is thus not good.
# In this case, increase relax to make smaller optimization steps, and also 
# maxiter to allow more iterations for convergence.

result = fit_peaks(ctx, backend=:qGN, relax=5)

## we print the result using
print_params(result.peak_results)

## and we plot the fit
plot_fit(ctx, result.peak_results)
savefig("fit_2.svg"); nothing #hide
# ![](fit_2.svg)

# You could also use the Optim backend, which uses a constrained L-BFGS-B search algorithm. 
# Therefore, in addition to using the prior uncertainties, 
# this backend leverages the lower and upper boundaries 
# declared earlier in `peaks_info`: 

result = fit_peaks(ctx, backend=:Optim)
print_params(result.peak_results)
plot_fit(ctx, result.peak_results)
savefig("fit_3.svg"); nothing #hide
# ![](fit_3.svg)

# ## Checking errors with bootstrapping
# The errors provided above come from the Hessian matrix at the optimal point.
# Those may not always be the good ones, if the loss function cannot be linearized close to the minimum.
#
# To check for parameter errors, one option is to use [bootstrapping](https://en.wikipedia.org/wiki/Bootstrapping_(statistics)).
# `bootstrap` allows to bootstrap your spectrum and refit the model on data subsamples.
# The interface is easy, similar to that of `fit_peaks` but the context is defined internally so 
# you don't even have to worry about that. Let's try it:

boot_params, boot_results = bootstrap(x, y_obs, estimated_error, peaks_info, nb_boot = 50, backend=:qGN, relax=5., maxiter=20)

# Here we only used 50 bootstraps for the sake of speed.
# In practice, you should use more bootstraps (e.g. 1000) to get a good estimate of the errors.
# The `bootstrap` function returns:
#   - a matrix of size (nb_params, nb_boot) with the fitted parameters (here `boot_params`);
#   - a peak_results Vector of Named Tuples (peak_type, parameters, areas), the values being tied to their errors thanks to Measurements.jl (here called `boot_results`).

# We can now print the bootstrapped results and compare the errors with those
# previously calculated from the Hessian:
print_params(boot_results)

# OK, actually for this example, we see that the errors from the boostrap analysis
# are close to those calculated from the Hessian matrix. Everything thus seems OK.
#
# Of course, a final quick visual inspection is always nice. This can be done by passing the 
# median of the matrix of bootstrapped parameters to the plot_fit function
plot_fit(ctx, boot_results)
savefig("fit_4.svg"); nothing #hide
# ![](fit_4.svg)

# ## Exploring errors on parameters using Turing.jl
#
# The same problem can be tackled using Turing.jl and the peak shape functions from Spectra as follow. 
# This offers another way to check that the estimated errors are good for instance.
#
# You will need to install Turing.jl, which is not a dependency of Spectra. 
# The code below runs well but it may not be fully optimized. It is just for the sack of example. 
# Also, it is commented for now because it is a bit too much for the documentation generation. 
# Run it on your own computer!
# If you have suggestions, do not hesitate!
# ```julia
# using Turing
#
# ## Define a Bayesian model with priors
# @model function bayesian_peaks(x, y)
#     ## Define priors based on peak_types
#     ## Example for a Gaussian peak
#     amplitude ~ truncated(Normal(10.016, 0.5), 0.0, Inf)
#     center ~ Normal(34.92, 0.5)
#     width ~ truncated(Normal(10.0, 0.5), 0.0, Inf)

#     μ = gaussian(x, [amplitude, center, width])
    
#     ## PEAK 2
#     amplitude2 ~ truncated(Normal(14.9, 0.5), 0.0, Inf)
#     center2 ~ Normal(55.0, 0.5)
#     width2 ~ truncated(Normal(3.0, 0.5), 0.0, Inf)
    
#     μ2 = lorentzian(x, [amplitude2, center2, width2])
    
#     ## PEAK 3
#     amplitude3 ~ truncated(Normal(25.5, 0.5), 0.0, Inf)
#     center3 ~ Normal(43.0, 10.0)
#     width3 ~ truncated(Normal(2.0, 0.5), 0.0, Inf)
#     lr ~ truncated(Normal(0.39, 0.03), 0.0, 1.0)
    
#     ## Calculate model prediction
#     μ3 = pseudovoigt(x, [amplitude3, center3, width3, lr])
    
#     ## Likelihood
#     σ ~ truncated(Normal(0.2, 0.03), 0.001, Inf)
#     y ~ MvNormal(μ + μ2 + μ3, σ^2 * I)
# end

# chain = sample(bayesian_peaks(x_fit, y_fit), NUTS(), 2000, nchains=3, progress=true)
# ```


