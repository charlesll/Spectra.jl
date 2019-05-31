# Pre-Processing

## Temperature and frequency corrections for Raman spectra

Raman spectra can be corrected from temperature and excitation line effects using this function.

```@docs
tlcorrection(data::Array{Float64},temp::Float64,wave::Float64;correction="long",normalisation="area",density=2210.0)
```

## Smoothing signals

Smoothing the signal is achieved with the smooth function.

```@docs
smooth(x,y;method="whittaker", window_length=5, polyorder = 2, Lambda = 10.0.^5, d=2, ese_y=1.0)
```

## Baseline subtraction

Baseline subtraction can be made with using the baseline function:

```@docs
baseline(x::Array{Float64},y::Array{Float64},roi::Array{Float64},basetype::AbstractString;polynomial_order=1, s = 1.0, lam = 10^5, p = 0.01, ratio = 0.01, niter = 10, p0_exp = [1.,1.,1.],p0_log =[1.,1.,1.])
```
Access to the gcvspline algorithm requires installation of the gcvspline Python package in Julia's
PyENV. This can be done like:

```julia-repl
 julia> using PyCall
 julia> pyimport_conda("pip", "pip")
 julia> run(`$(PyCall.python) -m pip install $(PACKAGES)`)
 ```

## Frequency shifts correction

In case your spectra are shifted from a reference value, Spectra offers several functions that allows you to correct it from this shift.

To correct a spectrum from a shift of P wavenumbers, you can simply call:

```@docs
xshift_direct(original_x::Array{Float64}, original_y::Array{Float64}, p::Float64)
```

Sometime, two signals from the same mineral show a shift in the X axis, while they share a common X axis. To correct from such thing, you can use the function:

```@docs
xshift_correction(full_x::Array{Float64}, full_shifted_y::Array{Float64}, ref_x::Array{Float64}, ref_y::Array{Float64},shifted_y::Array{Float64})
```

## array manipulation

For spectra recorded with decreasing frequencies, use the flipsp() function to
put them back with increasing frequencies (necessary for some algo)

```@docs
flipsp(spectra::Array{Float64})
```

You can also resample a signal at wanted x_new values with resample()

```@docs
resample(x::Array{Float64},y::Array{Float64},x_new::Array{Float64})
```
