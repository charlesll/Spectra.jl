# Measurements

Various functions to measure the parameters of peaks visible on spectra are available.

A new core function is introduced in v1.1: [`find_peaks`](@ref), which leverages the [Peaks.jl](https://www.juliapackages.com/p/peaks) package. This function allows you to, as its name says, find peaks in a signal and calculate various parameters, such as their heights, widths, centroids... The later calculation relies on the [`centroid`](@ref) function, which can be called independently.

Finally, earlier versions of Spectra had a [`peakmeas`](@ref) function to measure various parameters (height, width, centroid...). It works but please use [`find_peaks`](@ref) in your futur codes, it is much better. [`peakmeas`](@ref) will be removed in a future release of Spectra.

```@docs
find_peaks
centroid
peakmeas
```