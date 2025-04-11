# Measurements

Various functions to measure the parameters of peaks are available.

A new core function is introduced in v1.1: `find_peaks`, which leverages the [Peaks.jl](https://www.juliapackages.com/p/peaks) package:

```@docs
find_peaks
```

The centroid of a signal can be measured using the `centroid` function:

```@docs
centroid
```

Finally, I kept the old `peakmeas` function in the API for reference, but please use `find_peaks` in your futur codes, it is much better. `peakmeas` will be removed in a future release of Spectra.

```@docs
peakmeas
```