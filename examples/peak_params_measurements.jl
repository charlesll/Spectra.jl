# # Measuring peak parameters in a signal
#
# We first generate a dummy signal to play with:
#
# ## Signal generation
#
using Plots, Spectra, Random

## the x axis
x = collect(0:0.1:100)

## a fake signal: perfect y
y_perfect = gaussian(x, 1.0, 50.0, 6.0)

## we add noise: observed y
y = y_perfect + 0.05*randn(size(x,1))

plot(x, y, label="Noisy observed signal", xlabel="X", ylabel="Y")
savefig("ppm_1.svg"); nothing #hide
# ![](ppm_1.svg)

# ## Estimating parameters for one peak
#
# For one peak, estimations can be done using `peakmeas`:
#
height_y, hwhm_y, position_y, centroid_y, smoothed_y = peakmeas(x,vec(y),smoothing = "yes", method= "savgol",  y_smo_out = true)

println("Estimated peak height is $(height_y)")
println("Estimated peak hwhm is $(hwhm_y)")
println("Estimated peak position is $(position_y)")
println("Estimated peak centroid is $(centroid_y)")
# ## Estimating parameters on multiple peaks
#
# `peakmeas` is rather archaic. It does not work on multiple peaks. 
# To solve that, we introduced a new function, `find_peaks`, to do a better job:
#
result = find_peaks(x, y, smoothing=false, method="gcvspline", window_size = 20)
println("Peak positions: ", result.peak_positions)
println("Peak heights: ", result.peak_heights)
println("Peak hwhms: ", result.peak_hwhms)
println("Peak centroids: ", result.peak_centroids)
result.plot_peaks
savefig("ppm_2.svg"); nothing #hide
# ![](ppm_2.svg)

# Hum... we have multiple peaks? We need to tweak the parameters of the `find_peaks` function. 
# We can increase the `min_height` value to detect only peaks above a certain threshold. 
# We can also smooth the signal, it will help filter the noise.
# Other options are possible, see the documentation of `find_peaks`.

result = find_peaks(x, y, smoothing=true, method="gcvspline", window_size = 20, min_height=0.2)
println("Peak positions: ", result.peak_positions)
println("Peak heights: ", result.peak_heights)
println("Peak hwhms: ", result.peak_hwhms)
println("Peak centroids: ", result.peak_centroids)
result.plot_peaks
savefig("ppm_3.svg"); nothing #hide
# ![](ppm_3.svg)

# Good! It works now better.
#
# ## Peak centroid
#
# The centroid  of a peak or of a signal can also be measured using directly the `centroid` function. 
# It accepts x-y inputs, list of x-y spectra, or arrays of ys spectra associated to a vector of x values, see the documentation of `centroid`[@ref].
centroid2 = centroid(x,y)
println("Estimated peak centroid is $(centroid2)")
