# Smoothing signals

Smoothing the signal is achieved with the smooth function. Smoothing via the Whittaker, Savitsky-Golay and GCV cubic splines are available and provide good results: 

- The Whittaker implementation is a personal implementation, a convertion of the Matlab code of Eiler (2003). It can take in input a signal with x values equally spaced or not.
- The Savitsky-Golay implementation is from [SavitskyGolay.jl](https://github.com/lnacquaroli/SavitzkyGolay.jl) library. 
- The GCV cubic spline smoother behaves well, it is from the [DataInterpolations.jl library](https://github.com/SciML/DataInterpolations.jl)
- The Window-based smoother leverage the [DSP.jl library](https://github.com/JuliaDSP/DSP.jl).

The [`whittaker`](@ref) function can also be be used directly to have a fine control over the smoothing function, for instance by passing the weights `w` or changing `d` (also possible in `smooth`).

```@docs
smooth
whittaker
```