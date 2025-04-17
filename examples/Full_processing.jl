# # Processing data
#
# In this notebook, we reproduce the typical steps here for processing spectra
#    
# For this example, we create two Gaussian signals randomly sampled along two different X axis, with backgrounds
#
# ## Signal creation
using Spectra, Plots

## we create a fake signal with 
x_1 = rand(1000)*100
x_2 = rand(1000)*100

## create a signal that is the combination of two gaussian peaks plus a background
background_1 = 0.08 * x_1
background_2 = 0.03 * x_2

## some noise
noise_1 = 0.5*randn(1000)
noise_2 = 0.3*randn(1000)

y_1 = gaussian(x_1, 10.0, 40., 5.) .+ background_1 .+ noise_1
y_2 = gaussian(x_2, 20.0, 60., 9.) .+ background_2 .+ noise_2

## make a plot of our two spectra
scatter(x_1, y_1)
scatter!(x_2, y_2)

# ## First possible steps
# 
# We can do the following steps (not necessarily in this order):
#     
# - the X values are randomly sorted, we can solve that using [`flipsp`](@ref)
# - We may like to get our spectra on the same X axis for convenience using [`resample`](@ref)
# - Spikes are present: we can remove them using [`despiking`](@ref)
# - Backgrounds are present: we can remove them using [`baseline`](@ref)
# - Signals are noisy: we can smooth them using [`smooth`](@ref)
# - We can correct for temperature and laser wavelength using [`tlcorrection`](@ref)
# - We can normalize the spectra using [`normalise`](@ref)
# - We can fit the peaks using [`fit_peaks`](@ref)
#
# Let's do it!
#
# ## Sort X Axis
#
# We can sort the data by passing an array of spectra. After that we should have not problem plotting things with lines for instance!
spectrum_1 = flipsp([x_1 y_1])
spectrum_2 = flipsp([x_2 y_2])
plot(spectrum_1[:,1], spectrum_1[:, 2])
plot!(spectrum_2[:,1], spectrum_2[:, 2])
savefig("fp_1.svg"); nothing #hide
# ![](fp_1.svg)

# ## Resample to have everything on the same X axis
#
# Using [`resample`](@ref), we get everything on the same X axis. As we have two spectra with two different X axis, we simply provide them in a Vector like this:
x_new = collect(0.:0.5:100)
spectra_ = [[x_1 y_1], [x_2 y_2]]
spectra_same_x = resample(spectra_, x_new)
plot(x_new, spectra_same_x)
savefig("fp_2.svg"); nothing #hide
# ![](fp_2.svg)

# We see tiny problems with the interpolation, we can solve them using another method from DataInterpolations.jl, such as Linear:
spectra_same_x = resample(spectra_, x_new, method="LinearInterpolation")
plot(x_new, spectra_same_x)
savefig("fp_3.svg"); nothing #hide
# ![](fp_3.svg)

# ## Remove a background
# 
# Now we can remove a background using the [`baseline`](@ref) function. 
# Similarly to the other functions, you can pass x and y vectors or a x vectors and an array of y spectra. We will do that here:
ys_corrected, ys_baselines = baseline(x_new, spectra_same_x, method="arPLS")
p1 = plot(x_new, spectra_same_x)
plot!(x_new, ys_baselines)
savefig("fp_4.svg"); nothing #hide
# ![](fp_4.svg)

# ## Measure peak parameters
# 
# We can now measure the parameters of each peak, for instance their centroid:
println("Centroids are:")
println(centroid(x_new, ys_corrected))

# The peakmeas function to get access to a bunch of parameters: intensity, position, hwhm, centroïd
println("intensity, position, hwhm, centroïd are:")
println(peakmeas(x_new, ys_corrected[:,1]))

# If multiple peaks with a background are present, the best is to use the find_peaks() function
# It can even make a nice plot. Beware that you may have to tweak the window_size, min_width, min_height to detect only the 
# main peaks... See the documentation of [`find_peaks`](@ref) for details.
# for more precise usage, also see the Peaks.jl package.
result = find_peaks(x_new, ys_corrected[:,1], window_size=1, min_width=2., min_height=2.0, smoothing=true)
println("Peak positions: ", result.peak_positions)
println("Peak heights: ", result.peak_heights)
println("Peak HWHMs: ", result.peak_hwhms)
println("Peak centroids: ", result.peak_centroids)
result.plot_peaks

# ## Peak fitting!
#
# We can retrieve the parameters of the two peaks using peak fitting too! We do a quick fit with very loose 
# prior uncertainties on the peak parameters, and unrestrictive lower_bounds and upper_bounds.
# Those informations are provided in a vector of tuples, one tuple per peak:
# (type, initial_params, prior uncertainties, lower_bounds, upper_bounds)
peaks_info = [
        ## (type, initial_params, prior uncertainties, lower_bounds, upper_bounds)
        (:gaussian, [10.5, 50.0, 5.0], [100.0, 100.0, 100.0], [0., 0., 0.], [Inf, Inf, Inf]),
    ]

# we declare the context and fit the signals
ctx = prepare_context(x_new, peaks_info)
result_1 = fit_peaks(ctx, ys_corrected[:,1], backend=:Optim)
result_2 = fit_peaks(ctx, ys_corrected[:,2], backend=:Optim)

println("Parameters and fit for the first peak:")
print_params(result_1.peak_results, digits=3)
plot_fit(ctx, result_1.fitted_params, ys_corrected[:,1]; components=true)
savefig("fp_5.svg"); nothing #hide
# ![](fp_5.svg)

println("Parameters and fit for the second peak:")
print_params(result_2.peak_results, digits=3)
plot_fit(ctx, result_2.fitted_params, ys_corrected[:,2]; components=true)
savefig("fp_6.svg"); nothing #hide
# ![](fp_6.svg)

