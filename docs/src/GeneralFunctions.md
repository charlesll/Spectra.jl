# General Functions

## Peak shapes

The following functions are useful when generating peaks with various shapes. See the examples for using them during peak fitting for instance.

```@docs
gaussiennes(amplitude::Array{Float64},centre::Array{Float64},hwhm::Array{Float64},x::Array{Float64};style::String = "None")
lorentziennes(amplitude::Array{Float64},centre::Array{Float64},hwhm::Array{Float64},x::Array{Float64};style::String = "None")
pearson7(a1::Array{Float64},a2::Array{Float64},a3::Array{Float64},a4::Array{Float64},x::Array{Float64};style::String = "None")
pseudovoigts(amplitude::Array{Float64},centre::Array{Float64},hwhm::Array{Float64},lorentzian_fraction::Array{Float64},x::Array{Float64};style::String = "None")
normal_dist(nd_amplitudes::Array{Float64},nd_centres::Array{Float64},nd_sigmas::Array{Float64},x::Array{Float64})
```

## Peak measurement

```@docs
peakmeas(x::Array{Float64,1}, y::Array{Float64,1}; smoothing = "yes", method = "savgol", window_length=5, polyorder=2, ese_y=1., y_smo_out=false)
```

## Integration

Spectra.jl provides functions that allow one to integrate the area under a region of a spectrum, or to calculate the area under Gaussian, Lorentzian or other bands.

```@docs
trapz(x::Vector{Tx}, y::Vector{Ty}) where {Tx <: Number, Ty <: Number}
bandarea(Amplitude::Array{Float64},HWHM::Array{Float64}; peak_shape = "Gaussian", error_switch = "no", eseAmplitude::Array{Float64} = [0.0], eseHWHM::Array{Float64} = [0.0])
```

## Polynomial

```@docs
poly(p::Vector{Float64},x::Array{Float64})
```
