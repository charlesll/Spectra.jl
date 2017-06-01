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
peakmeas(x::Array{Float64},y::Array{Float64};smoothing = "yes", filter = :SavitzkyGolay, M=5,N=2,y_smo_out=false)
```

## Integration

Spectra.jl provides functions that allow one to integrate the area under a region of a spectrum, or to calculate the area under Gaussian, Lorentzian or other bands.

```@docs
trapz{Tx<:Number, Ty<:Number}(x::Vector{Tx}, y::Vector{Ty})
bandarea(Amplitude::Array{Float64},HWHM::Array{Float64}; peak_shape = "Gaussian", error_switch = "no", eseAmplitude::Array{Float64} = [0.0], eseHWHM::Array{Float64} = [0.0])
```

## Polynomials

```@docs
poly(p::Vector{Float64},x::Array{Float64})
polyfit(x::Array{Float64}, y::Array{Float64}, n::Int64)
```

## Splines

Not all the splines packages provide the same performances for data smoothing and interpolation. By experience, the Dierckx spline package ("Dspline" option in the baseline() function) provides a good starting point, but is not as usefull as other spline packages.

The csaps function of Matlab uses the SMOOTH Fortran library, and provides better smoothing capabilities for noisy data. Similarly, the GCVSPL Fortran package from Woltring (1986) also provides a very robust way to smooth and interpolate noisy data.

Starting from Spectra v0.3.4, the gcvspline Python module (https://github.com/charlesll/gcvspline) is used in the smooth and baseline function. It behaves exactly as the previous wrapping of GCVSPL.f in Julia, such that this should be transparent to users.

