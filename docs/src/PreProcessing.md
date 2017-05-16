# Pre-Processing

## Temperature and frequency corrections for Raman spectra

Raman spectra can be corrected from temperature and excitation line effects using this function.

```@docs
tlcorrection(data::Array{Float64},temp::Float64,wave::Float64;correction="long",normalisation="area",density=2210.0)
```

## Removing cristal or epoxy signals

Spectra.jl contains a function that helps removing the signal from crystals in the Raman spectra of glasses. Two spectra are needed: that of the mixed crystal+glass signals, and that of the pure cristal signals. Please note that it also can be used to remove signal from epoxy.

This function is still under test and experimental. Further details on the code will be provided soon. For now, only a short description is provided.

```@docs
	ctxremoval(liste,in_path,out_path,roi_all;input_properties=('\t',0),algorithm="FastICA",plot_intermediate_show = "no",plot_mixing_show = "yes",plot_final_show = "no",save_fig_switch = "yes", shutdown = 1300.,scaling=100.)
```

## Baseline subtraction

Baseline subtraction can be made with using the baseline function:

```@docs
baseline(x::Array{Float64},y::Array{Float64},roi::Array{Float64},basetype::AbstractString;p=1.0,SplOrder=3,roi_out="no")
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
	
## Smoothing

```@docs
SavitzkyGolayFilter{M,N}
```