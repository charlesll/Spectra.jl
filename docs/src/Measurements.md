# Measurements

## Introduction

Various functions to measure the parameters of peaks visible on spectra are available.

A new core function is introduced in v1.1: [`find_peaks`](@ref), which leverages the [Peaks.jl](https://www.juliapackages.com/p/peaks) package. This function allows you to, as its name says, find peaks in a signal and calculate various parameters, such as their heights, widths, centroids... The later calculation relies on the [`centroid`](@ref) function, which can be called independently.

Finally, earlier versions of Spectra had a [`peakmeas`](@ref) function to measure various parameters (height, width, centroid...). It works but please use [`find_peaks`](@ref) in your futur codes, it is much better. [`peakmeas`](@ref) will be removed in a future release of Spectra.

## Example

### Example signal

We first generate a dummy signal to play with:

```@example 1
using Plots, Spectra, Random

## the x axis
x = collect(0:0.1:100)

## a fake signal: perfect y
y_perfect = gaussian(x, 1.0, 50.0, 6.0)

## we add noise: observed y
y = y_perfect + 0.05*randn(size(x,1))

plot(x, y, label="Noisy observed signal", xlabel="X", ylabel="Y")
savefig("ppm_1.svg"); nothing #hide
```
![](ppm_1.svg)

### Estimating parameters for one peak

For one peak, estimations can be done using [`peakmeas`](@ref):

```@example 1
height_y, hwhm_y, position_y, centroid_y, smoothed_y = peakmeas(x,vec(y),smoothing = "yes", method= "savgol",  y_smo_out = true)

println("Estimated peak height is $(height_y)")
println("Estimated peak hwhm is $(hwhm_y)")
println("Estimated peak position is $(position_y)")
println("Estimated peak centroid is $(centroid_y)")
```

!!! warning

    `peakmeas` is basic. It does not work on multiple peaks. Prefer using [`find_peaks`](@ref).

### Find and measure multiple peaks

To find and measure parameters of a peak or multiple peaks in a signal, we can use [`find_peaks`](@ref):

```@example 1
result = find_peaks(x, y, smoothing=false, method="gcvspline", window_size = 20)
println("Peak positions: ", result.peak_positions)
println("Peak heights: ", result.peak_heights)
println("Peak hwhms: ", result.peak_hwhms)
println("Peak centroids: ", result.peak_centroids)
result.plot_peaks
savefig(result.plot_peaks, "ppm_2.svg"); nothing #hide
```
![](ppm_2.svg)

Hum... it finds too much peaks here. We need to tweak the parameters of the `find_peaks` function. 
We can increase the `min_height` value to detect only peaks above a certain threshold. 
We can also smooth the signal, it will help filter the noise.
Other options are possible, see the specific documentation of [`find_peaks`](@ref).

```@example 1
result = find_peaks(x, y, smoothing=true, method="gcvspline", window_size = 20, min_height=0.2)
println("Peak positions: ", result.peak_positions)
println("Peak heights: ", result.peak_heights)
println("Peak hwhms: ", result.peak_hwhms)
println("Peak centroids: ", result.peak_centroids)
result.plot_peaks
savefig(result.plot_peaks, "ppm_3.svg"); nothing #hide
```
![](ppm_3.svg)

Good! It works now better.

### Peak centroid

The centroid  of a peak or of a signal can also be measured using directly the [`centroid`](@ref) function. 
It accepts x-y inputs, list of x-y spectra, or arrays of ys spectra associated to a vector of x values, see the documentation of `centroid`[@ref].

```example 1
centroid2 = centroid(x,y)
println("Estimated peak centroid is $(centroid2)")
```

## Functions API

```@docs
find_peaks
centroid
peakmeas
```