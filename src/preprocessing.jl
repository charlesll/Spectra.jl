#############################################################################
#Copyright (c) 2016-2025 Charles Le Losq
#
#The MIT License (MIT)
#
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the #Software without restriction, including without limitation the rights to use, copy, #modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, #and to permit persons to whom the Software is furnished to do so, subject to the #following conditions:
#
#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, #INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
#############################################################################
"""
	xshift_direct(original_x::Array{Float64}, original_y::Array{Float64}, p::Float64)

Corrects a spectrum for a p shift in X. Used in xshift_correction.

# Inputs

	original_x: Array{Float64}
		x values
	original_y: Array{Float64}
		y values associated with x
	p: Array{Float64}
		value of how much x should be shifted

# Outputs

	original_x: Array{Float64}
		same as input
	corrected_y: Array{Float64}
		the y values corrected from the p shift in original_x
	p: Array{Float64}
		same as input.

"""
function xshift_direct(original_x::Array{Float64}, original_y::Array{Float64}, p::Float64)
    spl = Spline1D(original_x-p, original_y)
    corrected_y = spl(original_x)
    return original_x, corrected_y, p
end

"""
    correct_xshift(x::Vector{Float64}, y::Union{Vector{Float64}, Matrix{Float64}}, shift::Float64)
    correct_xshift(sps::Vector{<:Matrix}, shift::Float64)

Returns the signal(s) corrected from a given linear `shift` at the same `x` location as the input.

Signals can be provided as y (vector or an array of ys values) for a given x vector, or as a list of [x y] arrays of signals.

Depending on the arguments, it either returns a new vector or array of `y` at the position `x`, or a new list of corrected [x y] spectra.

This would be typically used to correct a linear shift in `x` on Raman spectra: for instance you measured the Si wafer peak at 522.1 cm-1 while you know it is at 520.7 cm-1. Therefore you will call this function to correct your spectra from this shift, without affecting the x values.

# Example

```julia
using Spectra

# for a vector y
x = [0., 1., 2., 3.]
y = 2*x
shift = -0.1

new_y = correct_xshift(x, y, shift)

# for an array of y
x2 = [0.5, 1.3, 2.0, 4.5]
y2 = [2*x 3*x 4*x]
new_y = correct_xshift(x, y2, shift)

# for a list of x-ys
old_spectra_list = [[x, y], [x2, y2]]
new_spectra_list = corrected(old_spectra_list, shift)
```
"""
function correct_xshift(x::Vector{Float64}, y::Vector{Float64}, shift::Float64)
    # method 1: y vector signal associated with x
    interp = AkimaInterpolation(
        y,
        x-shift;
        extrapolation_right=ExtrapolationType.Extension,
        extrapolation_left=ExtrapolationType.Extension,
    )
    return interp.(x)
end
function correct_xshift(x::Vector{Float64}, y::Matrix{Float64}, shift::Float64)
    # method 2: x vector and an array of y signals
    y_out = copy(y)
    for i in 1:size(y, 2)
        y_out[:, i] = correct_xshift(x, y[:, i], shift)
    end
    return y_out
end
function correct_xshift(sps::Vector{<:Matrix}, shift::Float64)
    # method 3: a list of spectra in 2 column arrays
    return [correct_xshift(sp[:, 1], sp[:, 2], shift) for sp in spectra]
end

"""

	flipsp(spectra::Union{Matrix{Float64}, Vector{<:Matrix}})

Returns the spectrum array or a list of spectra arrays sorted by their first column values (increasing x)

# Example

```julia

# create some unsorted signals
x = [0., 2.,5.,-1]
x2 = [0., 2.,5.,-1, 3., -10.]
y = 2*x
y2 = 2*x2

# the arrays of signals
signal_1 = [x y]
signal_2 = [x2 y2]

# a list of signal arrays
sps = [signal_1, signal_2]

# flip one signal array
flipsp(signal_1)

# you can also do this
flipsp([x y])

# flip the list of signal arrays
flipsp(sps)
```
"""
function flipsp(spectra::Matrix{Float64})
    x = spectra[:, 1]
    p = sortperm(x)
    return spectra[p, :]
end
function flipsp(spectra::Vector{<:Matrix})
    return [flipsp(sp) for sp in spectra]
end

"""
    resample(x::Vector{Float64}, y::Union{Vector{Float64}, Matrix{Float64}}, x_new::Vector{Float64}; method::String="AkimaInterpolation")
    resample(multiple_spectra::Vector{<:Matrix{Float64}}, x_new::Vector{Float64}; method::String="AkimaInterpolation")

Resample a signal or signals onto a new set of x-coordinates using interpolation.

# Arguments
- `x::Vector{Float64}`: The original x-coordinates corresponding to the input signal(s).
- `y::Union{Vector{Float64}, Matrix{Float64}}`: The input signal(s) to resample.
    - If `y` is a vector, it represents a single signal.
    - If `y` is a matrix, each column represents a separate signal.
- `multiple_spectra::Vector{<:Matrix{Float64}}`: A collection of spectra where each spectrum is a matrix with two columns:
    - First column: x-coordinates
    - Second column: y-values (signal intensities)
- `x_new::Vector{Float64}`: The new x-coordinates onto which the signal(s) will be resampled.
- `method::String="LinearInterpolation"`: The interpolation method to use. Options include:
    - methods available in the `DataInterpolations.jl` package: [https://docs.sciml.ai/DataInterpolations/stable/methods/](https://docs.sciml.ai/DataInterpolations/stable/methods/).

# Returns
- For single signal (vector) input: A vector of resampled values (`Vector{Float64}`).
- For multiple signals (matrix) input: A matrix where each column corresponds to the resampled values of the respective input column (`Matrix{Float64}`).
- For collection of spectra input: A matrix where each column contains the resampled y-values for the corresponding spectrum in the input collection.

# Errors
- Throws an error if an unsupported interpolation method is specified.
- Throws an error if the dimensions of `x` and `y` do not match.

# Methods
This function provides three methods to handle different input types:
1. Single signal: `resample(x::Vector{Float64}, y::Vector{Float64}, x_new::Vector{Float64})`
2. Multiple signals with common x-axis: `resample(x::Vector{Float64}, y::Matrix{Float64}, x_new::Vector{Float64})`
3. Collection of spectra: `resample(multiple_spectra::Vector{<:Matrix{Float64}}, x_new::Vector{Float64})`

# Notes
- Uses the `DataInterpolations.jl` package for interpolation.
- Extrapolation beyond the range of `x` is handled using the option `extrapolation_right=ExtrapolationType.Extension` and `extrapolation_left=ExtrapolationType.Extension`.
- For the collection of spectra method, each spectrum matrix must have exactly two columns: x-coordinates in the first column and y-values in the second column. However, something great: the different spectra can have different lengths!
- Automatically sorts the data

# Examples

## Example 1: resample a vector y or a matrix of ys (multiple spectra)
```@example
using Spectra, Plots

# signal creation
x = collect(0.:0.8:10.)
# create the signals with create_peaks()
peak_infos = [
    Dict(:type => :gaussian, :amplitude => 10.0, :center => 4.0, :hwhm => 0.6),
    Dict(:type => :gaussian, :amplitude => 5.0, :center => 6.0, :hwhm => 0.4),
]
ys, y = create_peaks(x, peak_infos)

# the new x axis
x_new = collect(0.:0.05:10.)

# resampling y as a vector
y2 = resample(x, y, x_new)

# resampling the ys array of the two peaks
y3 = resample(x, ys, x_new)

# plotting
p1 = scatter(x, y, label="Original data")
plot!(x_new, y2, label="Resampled y")
display(p1)

p2 = scatter(x, ys, label="Original peak data")
plot!(x_new, y3, label="Resampled peaks")
display(p2)
```

## Example 2: resampling a collection of spectra

```@example
x = collect(0.:0.8:10.)
y, ys = create_peaks(x, peak_infos)
x2 = collect(0.:0.8:10.)
y2, ys2 = create_peaks(x2, peak_infos)
x3 = collect(0.:0.8:10.)
y3, ys3 = create_peaks(x3, peak_infos)

spectra_ = [[x y], [x2 y2], [x3 y3]]
x_new = collect(0.:0.05:10.)
spectra_resampled = resample(spectra_, x_new)
p3 = scatter(x, [y, y2, y3], label="Original data")
plot!(x_new, spectra_resampled, label="Resampled data")
display(p3)
```

"""
function resample(
    x::Vector{Float64},
    y::Vector{Float64},
    x_new::Vector{Float64};
    method::String="LinearInterpolation",
)
    # this is the parent method: treat x-y Vectors
    if length(x) != length(y)
        throw(ArgumentError("x and y must have the same length"))
    end
    p = sortperm(x)# we automatically sort the data
    interp = getfield(DataInterpolations, Symbol(method))(
        y[p], x[p]; extrapolation = ExtrapolationType.Constant
    )

    return interp.(x_new)
end
function resample(
    x::Vector{Float64},
    y::Matrix{Float64},
    x_new::Vector{Float64};
    method::String="LinearInterpolation",
)
    # this method treats an array of y signals with a common X
    out = ones((size(x_new, 1), size(y, 2)))
    for i in 1:size(y, 2)
        out[:, i] = resample(x, vec(y[:, i]), x_new; method=method)
    end
    return out
end
function resample(
    multiple_spectra::Vector{<:Matrix{Float64}},
    x_new::Vector{Float64};
    method::String="LinearInterpolation",
)
    # this method treats a list of input spectra
    out = ones((size(x_new, 1), size(multiple_spectra, 1)))
    for (index, i) in enumerate(multiple_spectra)
        out[:, index] = resample(i[:, 1], i[:, 2], x_new; method=method)
    end
    return out
end

"""
    normalise(y::Union{Vector{Float64}, Matrix{Float64}}; x::Union{Vector{Float64}, Nothing}=nothing, method::String="intensity")

Normalise the y signal(s) using one of several methods.

# Arguments
- `y::Union{Vector{Float64}, Matrix{Float64}}`: The input signal(s) to normalize.
    - If `y` is a vector, it represents a single signal.
    - If `y` is a matrix, each column represents a separate signal.
- `x::Union{Vector{Float64}, Nothing}`: The x-coordinates corresponding to the y values (used only for `"area"` normalization). Default is `nothing`.
- `method::String`: The normalization method to use. Options are:
    - `"area"`: Normalize by the area under the curve (requires `x`).
    - `"intensity"`: Normalize by dividing by the maximum intensity.
    - `"minmax"`: Normalize to the range `[0, 1]`.

# Returns
- A normalized vector or matrix of signals (`Vector{Float64}` or `Matrix{Float64}`).

# Notes
- For `"area"` normalization, you must provide an x vector with the same length as each column of y.
- If using `"intensity"` or `"minmax"`, no x vector is required.

# Example

```julia
using Spectra, Plots

# Single signal normalization

x = collect(0.:0.1:10.)

# create a signal that is the combination of two gaussian peaks
y, ys = gaussiennes([10.,5.], [5.,6.], [1.,0.1], x)

# normalise the signal by area
y_norm = normalise(y; x=x, method="area")
plot(x, y_norm)

# Or normalise multiple signals, such as the two peaks created above (no need of the x axis in this case):

peaks_norm = normalise(ys, method="intensity")
plot(x, y_norm)
```

"""
function normalise(y::Vector{Float64}; x=nothing, method::String="intensity")
    # a function for Vector inputs
    if method == "area"
        if x == nothing
            throw(ArgumentError("Input array of x values for area normalisation"))
        end
        out = y ./ trapz(x, y)
    elseif method == "intensity"
        out = y ./ maximum(y)
    elseif method == "minmax"
        out = (y .- minimum(y)) ./ (maximum(y) - minimum(y))
    else
        throw(
            ArgumentError("Wrong method name, choose between area, intensity and minmax.")
        )
    end

    return out
end
function normalise(y::Matrix{Float64}; x=nothing, method::String="intensity")
    # a method for array inputs of Y values (several spectra, per column)
    out = copy(y)
    for i in 1:size(y)[2]
        out[:, i] = normalise(y[:, i]; x=x, method=method)
    end

    return out
end

"""
    despiking(x::Vector{Float64}, y::Vector{Float64}; neigh::Int=4, threshold::Int=3)
    despiking(x::Vector{Float64}, y::Matrix{Float64}; neigh::Int=4, threshold::Int=3)
    despiking(multiple_spectra::Vector{<:Matrix{Float64}}; neigh::Int=4, threshold::Int=3)

Remove spikes from signal(s) based on a threshold compared to a smoothed version.

This function smooths the spectra, calculates the residual error RMSE, and replaces points above 
threshold*RMSE with the average of non-spike neighboring points.

# Arguments
- `x::Vector{Float64}`: The x-coordinates of the signal.
- `y::Union{Vector{Float64}, Matrix{Float64}}`: The signal(s) to despike.
    - If `y` is a vector, it represents a single signal.
    - If `y` is a matrix, each column represents a separate signal.
- `multiple_spectra::Vector{<:Matrix{Float64}}`: A collection of spectra where each spectrum is a matrix with two columns:
    - First column: x-coordinates
    - Second column: y-values (signal intensities)
- `neigh::Int=4`: Number of points around each spike to consider for calculating the replacement value.
- `threshold::Int=3`: Multiplier of RMSE to identify spikes (points with residuals > threshold*RMSE).

# Returns
- For single signal input: A vector of despiked values (`Vector{Float64}`).
- For multiple signals input: A matrix where each column corresponds to the despiked values of the respective input column (`Matrix{Float64}`).
- For collection of spectra input: A vector of matrices, each containing the original x-coordinates and despiked y-values.

# Methods
This function provides three methods to handle different input types:
1. Single signal: `despiking(x::Vector{Float64}, y::Vector{Float64}; neigh::Int=4, threshold::Int=3)`
2. Multiple signals with common x-axis: `despiking(x::Vector{Float64}, y::Matrix{Float64}; neigh::Int=4, threshold::Int=3)`
3. Collection of spectra: `despiking(multiple_spectra::Vector{<:Matrix{Float64}}; neigh::Int=4, threshold::Int=3)`

# Notes
- The function uses the `smooth()` function with the "gcvspline" method to create a reference smoothed signal.
- Spikes are identified as points where the residual error exceeds threshold*RMSE.
- Spike values are replaced with the mean of neighboring non-spike points.

# Examples

## Example 1: Despiking a single signal
```@example
x = collect(0:0.1:10)
y = sin.(x) + 0.1*randn(length(x))
y[30] = 5.0 # Add a spike
y_clean = despiking(x, y)
```
## Example 2: Despiking multiple signals with common x-axis

```@example
x = collect(0:0.1:10)
y1 = sin.(x) + 0.1randn(length(x))
y2 = cos.(x) + 0.1randn(length(x))
y1[30] = 5.0 # Add a spike to first signal
y2[40] = -4.0 # Add a spike to second signal
y_matrix = hcat(y1, y2)
y_clean_matrix = despiking(x, y_matrix)
```

## Example 3: Despiking a collection of spectra
```@example
spectrum1 = hcat(collect(0:0.1:10), sin.(collect(0:0.1:10)) + 0.1randn(101))
spectrum1[30, 2] = 5.0 # Add a spike
spectrum2 = hcat(collect(0:0.1:8), cos.(collect(0:0.1:8)) + 0.1randn(81))
spectrum2[40, 2] = -4.0 # Add a spike
spectra_collection = [spectrum1, spectrum2]
clean_spectra = despiking(spectra_collection)
```

# Errors
- Throws an `ArgumentError` if `x` and `y` have different lengths.
- Throws an `ArgumentError` if `neigh` or `threshold` are not positive integers.
"""
function despiking(x::Vector{Float64}, y::Vector{Float64}; neigh::Int=4, threshold::Int=3)
    y_out = copy(y) # To avoid overwriting y
    # Check if x and y are of the same length
    if length(x) != length(y)
        throw(ArgumentError("x and y must have the same length"))
    end
    # Check if neigh is a positive integer
    if neigh < 0 || neigh != round(Int, neigh)
        throw(ArgumentError("neigh must be a positive integer"))
    end
    # Check if threshold is a positive integer
    if threshold < 0 || threshold != round(Int, threshold)
        throw(ArgumentError("threshold must be a positive integer"))
    end
    # smoothing
    y_smo = smooth(x, y; method="gcvspline")
    rmse_local = sqrt.((y - y_smo) .^ 2)
    rmse_mean = sqrt(mean((y - y_smo) .^ 2))

    # Identify spikes
    spikes = rmse_local .> threshold * rmse_mean

    for i in 1:length(y)
        if spikes[i] # If there is a spike in position i
            # Avoiding the boundaries
            low_i = max(1, i - neigh)
            high_i = min(length(y), i + 1 + neigh)

            w = low_i:high_i # Select 2m + 1 points around our spike
            w2 = w[.!spikes[w]] # Choose the ones which are not spikes
            y_out[i] = mean(y[w2]) # Average their values
        end
    end

    return y_out
end
function despiking(x::Vector{Float64}, y::Matrix{Float64}; neigh::Int=4, threshold::Int=3)
    # this method treats an array of y
    y_out = copy(y) # To avoid overwriting y
    # Check if x and y are of the same length
    if length(x) != size(y, 1)
        throw(ArgumentError("x and y must have the same length"))
    end
    # loop over y
    for i in 1:size(y, 2)
        y_out[:, i] = despiking(x, y[:, i]; neigh=neigh, threshold=threshold)
    end
    return y_out
end
function despiking(
    multiple_spectra::Vector{<:Matrix{Float64}}; neigh::Int=4, threshold::Int=3
)
    return [despiking(spectrum[:, 1], spectrum[:, 2]; neigh=neigh, threshold=threshold) for spectrum in multiple_spectra]
end

"""
    extract_signal(x::Vector{Float64}, y::Vector{Float64}, roi::Matrix{Float64}) -> Tuple{Vector{Float64}, Vector{Float64}, Vector{Int}}
    extract_signal(spectrum::Matrix{Float64}, roi::Matrix{Float64})
    extract_signal(multiple_spectra::Vector{<:Matrix{Float64}}, roi::Matrix{Float64})

Extract the `x`-`y` spectral values in specified regions of interest (ROI) and their indices.

You can pass associated x and y values, a single spectrum in the form of a [x y] matrix, or a list of [x y] matrices.

# Arguments
- `x::Vector{Float64}`: The x-axis values.
- `y::Vector{Float64}`: The corresponding y-axis values.
- `spectrum::Matrix{Float64}`: A matrix of size `n x 2`, where:
    - Column 1: x-axis values.
    - Column 2: y-axis values.
- `multiple_spectra::Vector{<:Matrix{Float64}}`: A collection of spectra, where each spectrum is a matrix with two columns:
    - First column: x-coordinates
    - Second column: y-values (signal intensities)
- `roi::Matrix{Float64}`: A matrix of size `n x 2`, where each row specifies a region of interest:
    - Column 1: Lower bounds of the ROI.
    - Column 2: Upper bounds of the ROI.

# Returns
- if calling `extract_signal(x, y, roi)`, it returns a tuple containing:
    1. `Vector{Float64}`: The x values within the ROI.
    2. `Vector{Float64}`: The y values within the ROI.
    3. `Vector{Int}`: The indices of the x and y values within the ROI.
- if calling `extract_signal(spectrum, roi)`, it returns a tuple containing:
    1. `Vector{Float64}`: The x values within the ROI.
    2. `Vector{Float64}`: The y values within the ROI.
    3. `Vector{Int}`: The indices of the x and y values within the ROI.
- if calling `extract_signal(multiple_spectra, roi)`, it returns a vector of tuples, where each tuple corresponds to a spectrum and contains:
    1. `Vector{Float64}`: The x values within the ROI.
    2. `Vector{Float64}`: The y values within the ROI.
    3. `Vector{Int}`: The indices of the x and y values within the ROI.

# Errors
- Throws an error if `x` and `y` have different lengths.
- Throws an error if `roi` is not a 2D matrix with exactly 2 columns.

# Notes
- The function sorts both `x`-`y` pairs and the ROI for safety before processing.
- Multiple ROIs are supported, and their results are concatenated.
"""
function extract_signal(x::Vector{Float64}, y::Vector{Float64}, roi::Matrix{Float64})
    # verify that x and y are the same length
    if length(x) != length(y)
        throw(ArgumentError("x and y must have the same length"))
    end
    # verify that roi is a 2D matrix
    if ndims(roi) != 2
        throw(ArgumentError("roi must be a 2D matrix"))
    end
    # verify that roi has 2 columns
    if size(roi, 2) != 2
        throw(ArgumentError("roi must have 2 columns"))
    end
    # verify that roi has at least 1 row
    if size(roi, 1) < 1
        throw(ArgumentError("roi must have at least 1 row"))
    end

    # sort x-y for safety
    p = sortperm(x)
    x_ = x[p]
    y_ = y[p]

    # sort roi for safety
    roi = sort(roi; dims=1)

    # find the indices of the x values that are in the roi
    interest_index::Vector{Int64} = findall(roi[1, 1] .<= x_ .<= roi[1, 2])
    if size(roi, 1) > 1
        for i in 2:size(roi, 1)
            interest_index = vcat(interest_index, findall(roi[i, 1] .<= x_ .<= roi[i, 2]))
        end
    end

    # we return the signal in interested regions and the index
    return x_[interest_index], y_[interest_index], interest_index
end
function extract_signal(spectrum::Matrix{Float64}, roi::Matrix{Float64})
    # method for an spectrum in an array [x y]
    return extract_signal(spectrum[:, 1], spectrum[:, 2], roi)
end
function extract_signal(multiple_spectra::Vector{<:Matrix{Float64}}, roi::Matrix{Float64})
    # method for a list of spectra
    extracted_signals = [extract_signal(spectrum, roi) for spectrum in multiple_spectra]
    return extracted_signals
end
function get_portion_interest(x::Vector{Float64}, y::Vector{Float64}, roi::Matrix{Float64})
    # method for an spectrum in an array [x y]
    println(
        "Old name function get_portion_interest() is now extract_signal(). Please use extract_signal() instead.",
    )
    return extract_signal(x, y, roi)
end

"""
    invcm_to_nm(shift_inv_cm::Vector{Float64}; laser_nm=532.0)

Converts Raman shifts in inverse centimeters (cm⁻¹) to absolute wavelengths in nanometers (nm), 
given `laser_nm`, the wavelength of the excitation laser in nanometers.

# Example
If using a 532 nm laser line, you will do:
```@example
x_inv_cm = collect(557:1.0:560) # unit = cm^-1
x_wavelength_nm = invcm_to_nm(x_inv_cm; laser_nm = 532.0)
```
"""
function invcm_to_nm(shift_inv_cm::Vector{Float64}; laser_nm::Float64=532.0)
    laser_inv_cm = 1.0 ./ (laser_nm .* 1.0e-7)
    x_inv_cm = laser_inv_cm .- shift_inv_cm
    x_nm = 1.0 ./ (x_inv_cm .* 1.0e-7)
    return x_nm
end

"""
    nm_to_invcm(x::Vector{Float64}; laser_nm = 532.0)

Converts absolute wavelengths in nanometers (nm) to Raman shifts in inverse centimeters (cm⁻¹), 
given `laser_nm`, the wavelength of the excitation laser in nanometers.

# Example
If using a 532 nm laser line, you will do:
```@example
x_wavelength_nm = collect(557:1.0:560) # unit = nm
x_inv_cm = nm_to_invcm(x_wavelength_nm; laser_nm = 532.0)
```
"""
function nm_to_invcm(x::Vector{Float64}; laser_nm::Float64=532.0)
    x_inv_cm = 1.0 ./ (x .* 1.0e-7)
    laser_inv_cm = 1.0 ./ (laser_nm*1.0e-7)
    shift_inv_cm = laser_inv_cm .- x_inv_cm
    return shift_inv_cm
end
