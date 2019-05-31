# Pre-Processing

## Temperature and frequency corrections for Raman spectra

Raman spectra can be corrected from temperature and excitation line effects using this function.

```@docs
tlcorrection
```

## Smoothing signals

Smoothing the signal is achieved with the smooth function.

```@docs
smooth
```

## Baseline subtraction

Baseline subtraction can be made with using the baseline function:

```@docs
baseline
```

## Frequency shifts correction

In case your spectra are shifted from a reference value, Spectra offers several functions that allows you to correct it from this shift.

To correct a spectrum from a shift of P wavenumbers, you can simply call:

```@docs
xshift_direct
```

Sometime, two signals from the same mineral show a shift in the X axis, while they share a common X axis. To correct from such thing, you can use the function:

```@docs
xshift_correction
```

## array manipulation

For spectra recorded with decreasing frequencies, use the flipsp() function to
put them back with increasing frequencies (necessary for some algo)

```@docs
flipsp
```

You can also resample a signal at wanted x_new values with resample()

```@docs
resample
```
