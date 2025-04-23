# Data Processing

## X axis correction and resampling

In case your spectra are shifted from a reference value, Spectra offers the [`correct_xshift`](@ref) function to correct it from this shift. 

For Raman spectroscopy, you can also sometime have a spectrometer returning an X axis in nm instead of cm``^{-1}``. In this case, convert the X axis from one unit to the other using either the [`nm_to_invcm`](@ref) or [`invcm_to_nm`](@ref) functions.

For spectra recorded with decreasing frequencies, use the [`flipsp`](@ref) function to
put them back with increasing frequencies (necessary for some algorithms).

Finally, if you want to have spectra sampled on the same X axis, or a spectrum sampled on a particular X axis, use the [`resample`](@ref) function.

```@docs
correct_xshift
nm_to_invcm
invcm_to_nm
flipsp
resample
```

## Spike removal

Spectra can be affected by spikes, sudden increase in signal intensity way above the sample signal. Those are also known as cosmic rays. The [`despiking`](@ref) function removes spikes / cosmic rays from spectra.

```@docs
despiking
```

## Temperature and frequency corrections for Raman spectra

For Raman spectroscopy, it is sometime useful to correct spectra from temperature and excitation line effects. The [`tlcorrection`](@ref) function does exactly this, and offers access to several mathematical functions available in the litterature.

```@docs
tlcorrection
```
## Signal normalisation and extraction

After some initial pre-processing, one may want to normalise the signal to its area or the maximum intensity. [`normalise`](@ref) can do that for you on a spectrum or even a vector of x-y spectra matrices.

You may be interested in specific portions of your signal. Extract them using [`extract_signal`](@ref).

```@docs
normalise
extract_signal
```

## Baseline subtraction

Baseline subtraction can be made with using the baseline function, which is the main API using other dedicated baseline functions. The documentation of all those functions is available below.

```@docs
baseline
```

The `baseline` function actually calls more dedicated functions, in particular for the automatic baselines. 

```@docs
als_baseline
arPLS_baseline
drPLS_baseline
```

The `rubberband` method does not work well for now but here is its docstrings for reference:

```@docs
rubberband_baseline
```
## Smoothing

Smoothing a signal or signals is achieved with the smooth function. Smoothing via the Whittaker, Savitsky-Golay and GCV cubic splines are available and provide good results: 

- The Whittaker implementation is a personal implementation, a convertion of the Matlab code of Eiler (2003). It can take in input a signal with x values equally spaced or not.
- The Savitsky-Golay implementation is from [SavitskyGolay.jl](https://github.com/lnacquaroli/SavitzkyGolay.jl) library. 
- The GCV cubic spline smoother behaves well, it is from the [DataInterpolations.jl library](https://github.com/SciML/DataInterpolations.jl)
- The Window-based smoother leverage the [DSP.jl library](https://github.com/JuliaDSP/DSP.jl).

The [`whittaker`](@ref) function can also be be used directly to have a fine control over the smoothing function, for instance by passing the weights `w` or changing `d` (also possible in `smooth`).

```@docs
smooth
whittaker
```
