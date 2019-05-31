# Peak fitting

## Model adjustment

Peak fitting is done with the JuMP framework (https://jump.readthedocs.org/en/latest/). Spectra.jl actually does not provide any peak fitting capacities by itself, but the combination of its functionality with JuMP helps making fitting procedure quite easy. An example is visible in the example section of Spectra.jl.

One goal of Spectra is to promote the use of global optimisation models, where peak parameters are actually calculated from variation in other parameters (chemistry, temperature, etc.), or are shared between several spectra. I will provide very soon an example of such an approach. It can be implemented in a few lines of code with combining Spectra and JuMP, and has the advantage of greatly reducing the errors of the fits.

## error calculation with bootstrapping

Error calculation can be done with using bootstrapping. Spectra provides a function that allows generating K new datasetes, by resampling the existing dataset in a non-parametric or parametric way.

```@docs
bootsample(x::Array{Float64}, y::Array{Float64}; boottype::String = "np", ese::Array{Float64} = [0.0])
bootperf(params_boot::Array{Float64}; plotting::String = "True", parameter::Int64 = 1, feature::Int64 = 1, histogram_step::Int64 = 100, savefigures::String = "False", save_bootrecord::String = "Boot_record.pdf", save_histogram::String = "Boot_histogram.pdf")
```

For further details, see the following references

Efron, B. 1979. “Bootstrap Methods: Another Look at the Jackknife.” The Annals of Statistics 7 (1): 1–26.

Efron, Bradley. 1981. “Nonparametric Estimates of Standard Error: The Jackknife, the Bootstrap and Other Methods.” Biometrika 68 (3): 589–99. doi:10.1093/biomet/68.3.589.

Efron, B., and Tibshirani, R. 1994. An Introduction to the Bootstrap. CRC press.
