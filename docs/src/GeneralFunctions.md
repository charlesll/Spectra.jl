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
peakhw(x::Array{Float64},y::Array{Float64};M=5,N=2,y_smo_out=false)
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

This GCVSPL spline package is called directly by Julia (through a ccall()) in the baseline function, with the options of a cubic spline with least-square data fitting. The smoothing is done with scaling the variances of the data points (VAR variable in the GCVSPL.f package) that is provided to the GCVSPL.f program.

Now, while baseline() should be well suited for most users needs, it uses cubic splines that are not always the best answers to some problems. For instance, quadratic splines may be more robust in some cases. You can change that by providing the spline order to baseline() as SplOrder = 2 for instance.

In case you want to have even more control on GCVSPL.f, and use its internal tricks and tweeks, the following lines will provide you the documentation of the two functions allowing you to calculate the spline coefficients and to evaluate the spline values at specific x entries.

```@docs
gcvspl_julia(x::Array{Float64,1},y::Array{Float64,1},ese::Array{Float64,1},SmoothSpline::Float64;SplineOrder::Int32 = Int32(2),SplineMode::Int32 = Int32(3))

splderivative_julia(xfull::Array{Float64},xparse::Array{Float64},cparse::Array{Float64};SplineOrder::Int32 = Int32(2), L::Int32 = Int32(1), IDER::Int32 = Int32(0))

```
