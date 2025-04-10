# Peak fitting

## Model adjustment

Spectra provides a peak fitting function to perform simple fits with Gaussian, Loretnzian and Pseudovoigt peaks, using Optim.jl.

```@docs
fit_peaks
```

## Interaction with JuMP.jl

For more complex fits, combining the abilities of Spectra with JuMP helps making fitting procedure quite easy. An example is visible in the example section of Spectra.

One goal of Spectra is to promote the use of global optimisation models, where peak parameters are actually calculated from variation in other parameters (chemistry, temperature, etc.), or are shared between several spectra. I will provide very soon an example of such an approach. It can be implemented in a few lines of code with combining Spectra and JuMP, and has the advantage of greatly reducing the errors of the fits.

## Error calculation with bootstrapping

Error calculation can be done with using bootstrapping. Spectra provides a function that allows generating K new datasetes, by resampling the existing dataset in a non-parametric or parametric way.

```@docs
bootsample
bootperf
```

