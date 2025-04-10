# Pre-Processing

## Temperature and frequency corrections for Raman spectra

Raman spectra can be corrected from temperature and excitation line effects using this function.

```@docs
tlcorrection
```

## Frequency shifts correction

In case your spectra are shifted from a reference value, Spectra offers several functions that allows you to correct it from this shift.

```@docs
xshift_direct
correct_xshift
```

## Array manipulation

For spectra recorded with decreasing frequencies, use the flipsp() function to
put them back with increasing frequencies (necessary for some algo)

```@docs
flipsp
```

## Resampling

```@docs
resample
```

## Spike removal

We have a function to remove spikes from spectra.

```@docs
despiking
```

## Signal normalisation

```@docs
normalise
```
## Signal extraction

```@docs
extract_signal
```