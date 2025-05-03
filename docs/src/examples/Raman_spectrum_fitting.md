The source files for all examples can be found in [/examples](https://github.com/charlesll/Spectra.jl/tree/master/examples/).
```@meta
EditURL = "../../../examples/Raman_spectrum_fitting.jl"
```

# Peak fitting the Raman spectrum of a glass

Code by Charles Le Losq, Created 7 April 2015 for Python, Modified 30 Sept. 2016 for Julia, updated February 2019, then April 2025.

Last modified: April 2025 for release of Spectra v2.0.0.

Glass materials usually present broad bands due to their inherent structural disorder.
Due to this, peak-fitting Raman spectra of glasses is informative but difficult, and requires
a careful balance between model freedom and constraints. Using boundaries for parameters can help,
as well as setting priors on model parameters and working in a Bayesian framework. Spectra provides a function
to do such things, `fit_peaks`.

Here, we will fit a Raman spectrum of a glass with Julia, using the `fit_peaks` function now available in Spectra.
We will leverage the quasi-Newton method, with a loss function that takes into account data and prior model errors.
In other terms, we are setting ourselves in a Bayesian framework, with priors on the model parameters.

You could solve the following problem in a pure Bayesian approach, using e.g. HMC algorithm in Turing.jl.
However, it takes some time... The quasi-Newton method is clean, the Julia code is the direct transcription of the mathematical formulas available in Tarantola ([2005](https://epubs.siam.org/doi/book/10.1137/1.9780898717921)), chapter 3.
This method is fast, accurate, but sometimes a bit instable.
You can also use the Interior Point Newton algorithm from Optim.jl, by setting the backend to `:Optim`. It is good, stable but slightly slower. It also benefits from boundaries!

## Context

In this example, we fit the 850-1300 cm$^{-1}$ portion of a Raman spectrum of a lithium tetrasilicate glass Li$_2$Si$_4$O$_9$,
the name will be abbreviated LS4 in the following.

For further references for fitting Raman spectra of glasses, please see for instance:

- Virgo et al., 1980, Science 208, p 1371-1373;
- Mysen et al., 1982, American Mineralogist 67, p 686-695;
- McMillan, 1984, American Mineralogist 69, p 622-644; Mysen, 1990, American Mineralogist 75, p 120-134;
- Le Losq et al., 2014, Geochimica et Cosmochimica Acta 126, p 495-517
- Le Losq et al., 2015, Progress in Earth and Planetary Sciences 2:22.

`fit_peaks` is meant to fit a single spectrum. You could do a loop to fit several spectra, but if you want to define global models, I invite you to use the [JuMP framework](https://jump.dev/) or, if you fancy Bayesian MCMC methods, [Turing.jl](https://turinglang.org/index.html)

## Importing libraries and data

First, we import the libraries for doing various things:

````@example Raman_spectrum_fitting
using Spectra ## our Spectra library
using Statistics ## to have access to core functions like mean() or std()
using DelimitedFiles ## to import the data

# Plotting libraries
using Plots
gr()
# LaTeX for superscripts/subscripts in labels, captions...
using LaTeXStrings
````

## Importing data

For that, we use readdlm from DelimitedFiles, and we skip the header and footer lines.

````@example Raman_spectrum_fitting
data = readdlm(joinpath(@__DIR__, "data/LS4.txt"), '\t', Float64)
skip_header = 23
skip_footer = 121
inputsp = zeros(size(data)[1]-skip_header-skip_footer,2)
j = 1
for i = skip_header+1:size(data)[1]-skip_footer
    inputsp[j,1] = Float64(data[i,1])
    inputsp[j,2] = Float64(data[i,2])
    global j += 1
end
````

We will now make the following pre-processing: we correct the data from temperature and excitation line effects using the `tlcorrection` function.
It is not always necessary at frequencies > 500 cm-1,
but this is just for the sack of example in the present case

````@example Raman_spectrum_fitting
x, y, ese_y = tlcorrection(inputsp, 23.0, 490.0)
````

and now we create a new plot for showing the spectrum

````@example Raman_spectrum_fitting
plot(x, y,
xlabel=L"Raman shift, cm$^{-1}$",
ylabel="Normalized intensity, a. u.",
title="Figure 1: the spectrum of interest")
savefig("rsf_1.svg"); nothing #hide
````

![](rsf_1.svg)

So we are looking at the 500-1300 cm$^{-1}$ portion of the Raman spectrum of the glass.
We see a peak near 800 cm$^{-1}$, and two others near 950 and 1085 cm$^{-1}$.
We will be interested in fitting the 870-1300 cm$^{-1}$ portion of this spectrum, which can be assigned to
the various symmetric and assymetric stretching vibrations of Si-O bonds in the SiO$_2$ tetrahedra
present in the glass network (see the above cited litterature for details).

## Baseline Removal

First thing we notice in Fig. 1, we have to remove a baseline because this spectrum is shifted from 0 by some "background" scattering. This quite typical in Raman spectra of glasses.
Several ways exist to do so. We're going to the simplest thing: a polynomial fitting the signal base around 870 and 1300 cm$^{-1}$.
Other reasonnable solutions include a linear function, and a constant function. The two latter can be fitted between 1300 and 1350 cm$^{-1}$,
but we will need to add another peak around 800 cm$^{-1}$. For now, the example is done with fitting the 870 cm$^{-1}$ portion of spectra,
as this usually results in more robust final results.

First, we define the regions of interest roi where we think the baseline is.
Then, we call the `baseline` function to define the baseline and subtract it from y.

````@example Raman_spectrum_fitting
roi = [860.0 870.0; 1300.0 1400.0]
y_corr, y_bas = baseline(x, y, roi=roi, method="polynomial", polynomial_order=2)

# To visualize this, we create a plot showing the baseline:
plot(x, [y y_corr y_bas],
xlabel=L"Raman shift, cm$^{-1}$",
ylabel="Normalized intensity, a. u.",
title="Figure 2: the fit of the background")
savefig("rsf_2.svg"); nothing #hide
````

![](rsf_2.svg)

## Signal extraction and normalisation

Now we will get the portion of the spectrum we want to fit using `extract_signal`. Then we normalise the signal using `normalise`.
We will calculate the errors based on the comparison between the signal and its "smoothed version", provided by `smooth`.

````@example Raman_spectrum_fitting
# First we extract the signal we want to fit:
x_fit, y_fit, _ = extract_signal(x, y_corr, [860. 1300.])

# We normalise y_fit so that its area is 1
y_fit = normalise(y_fit; x=x_fit, method="area")

# We smoothed the signal and get an estimate of the errors using it.
y_fit_perfect = smooth(x_fit, y_fit, method="whittaker", lambda=1e2)
ese_y_fit = sqrt.(mean((y_fit_perfect-y_fit).^2)) * ones(size(y_fit_perfect))

# Let's have a look at the signal in the fitting region
plot(x_fit, y_fit, label="Signal to fit")
plot!(x_fit, y_fit_perfect, label="smoothed",
    xlabel=L"Raman shift, cm$^{-1}$",
    ylabel="Normalized intensity, a. u.",
    title="Figure 3: signal to peak fit")
savefig("rsf_3.svg"); nothing #hide
````

![](rsf_3.svg)

## Fitting the spectrum

We do the fit using the `fit_peaks` function with the quasi-Newton algorithm. From the litterature, we have five peaks (see Le Losq et al. 2014, 2015 and references therein).
Compared to earlier studies, we also now know that the main one near 1080 cm$^-1$ may be actually a pseudovoigt peak.

Here we will place a strong prior on the intensity of the peak near 1090 cm$^-1$. It seems even a bit unrealistic but the influence of the prior loss compared
to the data loss is somehow small, so if you want to add tight constraints you sometime need to place low uncertainties on the prior values of the parameters you want to constrain.
Here we assume that the central band should be fit by a peak with a strong intensity. We use a 1 \% prior uncertainty on a prior value that actually is 10 \% lower than the maximum intensity:

````@example Raman_spectrum_fitting
max_y = maximum(y_fit)
prior_main_peak = max_y-0.1*max_y
prior_main_peak_uncertainty = 0.01*prior_main_peak
````

Let's implement that in a `peaks_info` vector of peak parameters, uncertainties and boundaries.

Then we declare the context and have a look at the prior mode. If necessary we re-adjust it.
Keep in mind that it should be fairly close to the solution, as the algorithms (here the quasi-Newton method)
we use are local search algorithms.

The `peaks_info` vector is a list of tuples, each tuple containing the following information:
- the type of peak (e.g. :gaussian, :lorentzian, :pseudovoigt)
- the parameters of the peak (e.g. amplitude, position, width, shape)
- the uncertainties on the parameters (e.g. amplitude, position, width, shape)
- the lower bounds for the parameters (e.g. amplitude, position, width, shape)
- the upper bounds for the parameters (e.g. amplitude, position, width, shape)

````@example Raman_spectrum_fitting
peaks_info = [
        (:gaussian,    [0.002, 950, 27],       [0.0005, 5.0,3.0], [0.0, 0.0, 0.0], [Inf, Inf, 60.0]),
        (:gaussian,    [0.0044, 1044, 40],      [0.0005, 5.0, 3.0], [0.0, 0.0, 0.0], [Inf, Inf, 60.0]),
        (:pseudovoigt,   [prior_main_peak, 1086, 30., 0.8], [prior_main_peak_uncertainty, 5.0, 3.0, 0.02], [0.0, 0.0, 0.0, 0.0], [Inf, Inf, 60.0, 1.0]),
        (:gaussian,    [0.0028, 1150, 45],      [0.0005, 5.0, 3.0], [0.0, 0.0, 0.0], [Inf, Inf, 60.0]),
        (:gaussian,    [0.0009, 1185, 30],      [0.0001, 5.0, 3.0], [0.0, 0.0, 0.0], [Inf, Inf, 60.0]),
    ]

# we declare the context of the fit
ctx = prepare_context(x_fit, y_fit, peaks_info, ese_y_fit)

# We plot the prior model
p = plot_fit(ctx, title="Prior model")
savefig("rsf_4.svg"); nothing #hide
````

![](rsf_4.svg)

````@example Raman_spectrum_fitting
# doing the fit
result = fit_peaks(ctx, backend=:Optim, relax=5, maxiter=100)

# we print the result using
print_params(result.peak_results)

# and we plot the fit
plot_fit(ctx, result.peak_results)
savefig("rsf_5.svg"); nothing #hide
````

![](rsf_5.svg)

## Checking errors with bootstrapping

The error bars above appear fairly small... We check them using boostrapping.
Again, we use a unrealistically low number of bootstrapped samples here because this code runs during documentation generation,
but in reality you would like to increase `nb_boot` to something like 1000.
I also fidled with `maxiter` to set it to a small value while still seeing convergence, such that we don't spend too much time in the quasi-Newton algorithm.

````@example Raman_spectrum_fitting
boot_samples, boot_results = bootstrap(x_fit, y_fit, ese_y_fit, peaks_info, nb_boot = 50, backend=:qGN, relax=5., maxiter=20)
````

We can now print the bootstrapped results and compare the errors with those
previously calculated from the Hessian:

````@example Raman_spectrum_fitting
print_params(boot_results)
````

OK, actually for this example, we see that the errors from the boostrap analysis
are close to those calculated from the Hessian matrix. Everything thus seems OK.

Of course, a final quick visual inspection is always nice. This can be done by passing `boot_results`
to the `plot_fit` function

````@example Raman_spectrum_fitting
plot_fit(ctx, boot_results)
savefig("rsf_6.svg"); nothing #hide
````

![](rsf_6.svg)

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

