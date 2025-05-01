# Data Processing

## Introduction

Spectra allows you to perform several processing steps on x-y spectral data. Below we will showcase short examples, and then you will find the documentation of the variosu functions you may want to use!

As a starting point and for the sack of example, we create two synthetic signals to play with. They will be Gaussian signals randomly sampled along two different X axis, with noise and increasing backgrounds. One of them will also have a strong spike!

```@example 1
# Signal creation
using Spectra, Plots

# we create a fake signal with 
x_1 = rand(1000)*100
x_2 = rand(1000)*100

# create a signal that is the combination of two gaussian peaks plus a background
background_1 = 0.08 * x_1
background_2 = 0.03 * x_2

# some noise
noise_1 = 0.5*randn(1000)
noise_2 = 0.3*randn(1000)

# the full toy signals
y_1 = gaussian(x_1, 10.0, 40., 5.) .+ background_1 .+ noise_1
y_2 = gaussian(x_2, 20.0, 60., 9.) .+ background_2 .+ noise_2

# one of them will have a spike!
y_1[600] = 250.0 

# make a plot of our two spectra
scatter(x_1, y_1)
scatter!(x_2, y_2)
savefig("fp_1.svg"); nothing #hide
```
![](fp_1.svg)

We can do the following steps (not necessarily in this order):

- [`correct_xshift`](@ref) allows correcting X-axis shifts of your spectra from a reference value (e.g. silicon wafer reference in Raman spectroscopy). 
- [`nm_to_invcm`](@ref) or [`invcm_to_nm`](@ref) convert the X axis between nanometers (nm) and wavenumbers (cm``^{-1}``).
- [`flipsp`](@ref) sort the X-axis (this is necessary for some algorithms).
- [`resample`](@ref) allows getting our spectra on the same X axis for convenience. 
- [`despiking`](@ref) remove spikes in the signal. 
- [`baseline`](@ref) allows removing the background.
- [`smooth`](@ref) allows smoothing signals.
- [`tlcorrection`](@ref) corrects Raman spectra for temperature and laser wavelength effects.
- [`normalise`](@ref) allows normalising the spectra. 
- [`extract_signal`](@ref) can extract specific portions of a signal.

Thanks to Julia's multiple dispatch, those functions support different types of inputs. Of course you will receive different outputs, see the individual documentation of each function for further details. This is quite convenient as this avoid you to write your own loops to process many spectra at once.

Let's now use some of those functions below on the signal generated above.

## Sort X Axis

We can sort the data by passing an array of spectra to [`flipsp`](@ref). After that we should have not problem plotting things with lines for instance!

```@example 1
spectrum_1 = flipsp([x_1 y_1])
spectrum_2 = flipsp([x_2 y_2])
plot(spectrum_1[:,1], spectrum_1[:, 2])
plot!(spectrum_2[:,1], spectrum_2[:, 2])
savefig("fp_2.svg"); nothing #hide
```
![](fp_2.svg)

## Remove spikes

OK, the plot above reveals a strong spike in one of the signals. We will treat actually both signals with
[`despiking`](@ref) to remove possible spikes from the signals, using the default parameters.  
In summary, with the default settings, [`despiking`](@ref) checks if any points in a spectrum differ by more than 3 sigma from the mean value of the 4 neighboring points.  
You can change the default values to adjust the threshold (for more or less than 3-sigma), or to modify the number of neighboring points considered.

```@example 1
y_1 = despiking(x_1, y_1)
y_2 = despiking(x_2, y_2)
```

!!! tip

    You could also call `despiking` on the collection of spectra as
    ```julia
    collection_spectra = [[x_1 y_1], [x_2 y_2]] 
    ys = despiking(collection_spectra)
    ```

## Resample spectra

Using [`resample`](@ref), we can resample a spectrum or spectra on a user-defined X axis by calling

```julia
x_new = collect(0.:0.5:100)
y_new = resample(x, y, x_new)
```

By default, [`resample`](@ref) uses a linear interpolation method from the `DataInterpolations.jl` package, but you can specify other methods available at [https://docs.sciml.ai/DataInterpolations/stable/methods/](https://docs.sciml.ai/DataInterpolations/stable/methods/).

If you have multiple spectra, it is here very interesting to provide a collection of those spectra because you will then receive
an array of spectra in output, all sampled on the same X axis.

Continuing on the example shown above, we can do:
```@example 1
x_new = collect(0.:0.5:100)
spectra_ = [[x_1 y_1], [x_2 y_2]]
spectra_same_x = resample(spectra_, x_new)
plot(x_new, spectra_same_x)
savefig("fp_3.svg"); nothing #hide
```
![](fp_3.svg)

## Baseline subtraction

Baseline subtraction is performed using [`baseline`](@ref), which serves as the main API and wraps several dedicated baseline correction algorithms. Similarly to the other functions, you can pass x and y vectors or a x vectors and an array of y spectra. 

Continuing with our example, we will do here:
```@example 1
ys_corrected, ys_baselines = baseline(x_new, spectra_same_x, method="arPLS")
p1 = plot(x_new, spectra_same_x)
plot!(x_new, ys_baselines, labels=["background 1" "background 2"])
savefig("fp_4.svg"); nothing #hide
```
![](fp_4.svg)

Other methods are available, see the Tutorials and [`baseline`](@ref) function documentation  for further details!

## Smoothing

Spectra smoothing can be achieved with the [`smooth`](@ref) function, which supports several algorithms: 

- **Whittaker smoother**: Custom Julia implementation based on the Matlab code of Eiler (2003). It supports both equally and unequally spaced X values.
- **Savitzky-Golay Smoother**: Provided by the [SavitskyGolay.jl](https://github.com/lnacquaroli/SavitzkyGolay.jl) library. 
- **GCV cubic spline smoother:** From the [DataInterpolations.jl](https://github.com/SciML/DataInterpolations.jl) library.
- **Window-based smoothers**: leverage the [DSP.jl](https://github.com/JuliaDSP/DSP.jl) library.

For fine control over smoothing parameters, you can use the [`whittaker`](@ref) function directly, allowing you to change weights `w` or the smoothing order `d` (also possible in `smooth`).

Continuing with our example, we will pass the matrix of baseline corrected signals to `smooth` like:
```@example 1
smoothed_y = smooth(x_new, spectra_same_x; method="gcvspline")

p1 = plot(x_new, spectra_same_x)
plot!(x_new, smoothed_y, labels=["smoothed 1" "smoothed 2"])
savefig("fp_5.svg"); nothing #hide
```
![](fp_5.svg)

Other methods are available, see the Tutorials and [`smooth`](@ref) function documentation  for further details!

## Signal normalisation

Using [`normalise`](@ref), you can normalise signals to their maximum intensity (`method="intensity"`), the area under the curve (`method="area"`) or to their minimum and maximum values (minimum will be set to 0, maximum to 1) (`method="minmax"`).

For instance, continuing with our example, we can do:

```@example 1
normalised_ys = normalise(spectra_same_x, method="intensity")
p1 = plot(x_new, normalised_ys)
savefig("fp_6.svg"); nothing #hide
```
![](fp_6.svg)

If you want to normalize the signals by their areas, you have to pass `x` values too, like:

```julia
normalised_y = normalise(y_matrix, x, method="intensity")
```

## Signal extraction

Extract signals in specific regions of interest using [`extract_signal`](@ref). You can pass associated x and y values, a single spectrum in the form of a [x y] matrix, or a list of [x y] matrices.

For instance, for a single signal in which we want the values between 40 and 60, we would write:

```julia
roi = [[40. 60.]]
extracted_x, extracted_y, indices = extract_signal(x, y, roi)
```

You can also extract the signals in different portions by using a matrix for the regions of interest. For instance, to extract signals between 20 and 40, and 60 and 80, we can do:
```julia
roi = [[20. 40.]; [60. 80.]]
extracted_x, extracted_y, indices = extract_signal(x, y, roi)
```

## Functions API

```@docs
correct_xshift
nm_to_invcm
invcm_to_nm
flipsp
resample
despiking
tlcorrection
normalise
extract_signal
baseline
als_baseline
arPLS_baseline
drPLS_baseline
rubberband_baseline
smooth
whittaker
```
