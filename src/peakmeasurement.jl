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
# functions.jl contains several mathematic functions
#
#
#############################################################################

"""

    peakmeas(x::Array{Float64,1}, y::Array{Float64,1}; smoothing::Bool = true, method::String = "savgol", window_length::Int=5, polyorder::Int=2, ese_y::Float64=1., y_smo_out::Bool=false)

Perform measurements of the position, width, intensity and centroid of a dominant peak in a provided `x`-`y` signal.

It smooths the signal with a Savitzky-Golay filter prior to measuring the peak position, 
width and intensity. It is advised to check that the M and N values of the Savitzky-Golay 
filter are adequate for your problem before trusting the results from peakmeas. For that, set `y_smo_out=true`.

Half-width at half-maximum are calculated as the width of the peak at half its maximum intensity. This calculation is not affected by any asumption of peak symmetry (no fitting is done).

Centroïd is calculated as sum(y./sum(y).*x).

!!! warning

    This function may not stay in future versions. Consider using `find_peaks` instead.

# Inputs
-`x::Array{Float64}`: the x values
-`y::Array{Float64}`: the y values

# Options
-`smoothing::String`: triggers the smoothing of the spectrum if set to yes (default value);
-`filter::Symbol`: the filter that will be used. See the smooth function documentation;
-`M::Int`: M parameter for smoothing y with a Savitzky-Golay filter. See smooth function documentation. Default = 5.
-`N::Int`: N parameter for smoothing y with a Savitzky-Golay filter. See smooth function documentation. Default = 2. 
-`y_smo_out::bool`: Outputs the smoothed signal.

# Outputs
-`intensity::Float64`: peak intensity
-`position::Float64`: peak position
-`hwhm: Float6:: peak` half-width at half-maximum
-`centroïd::Float64`: peak centroid
"""
function peakmeas(
    x::Vector{Float64},
    y::Vector{Float64};
    smoothing="yes",
    method="savgol",
    window_length=5,
    polyorder=2,
    ese_y=1.0,
    y_smo_out=false,
)
    ### PRELIMINARY CHECK: INCREASING SIGNAL
    sp = flipsp([x y])
    x = deepcopy(sp[:, 1])
    y = deepcopy(sp[:, 2])

    if smoothing == "yes"
        y_smo = smooth(
            x, y; method=method, window_length=window_length, polyorder=polyorder
        )
    else
        y_smo = collect(y)
    end

    x_maximum = x[y_smo .== maximum(y_smo)]

    if size(x_maximum, 1) >= 2
        error(
            "Something went wrong, is there two peaks with identical intensities in the signal? This function treats signal with one main peak...",
        )
    end

    x_1 = x[x .< x_maximum];
    x_2 = x[x .>= x_maximum];
    y_first_portion = y_smo[x .< x_maximum];
    y_second_portion = y_smo[x .>= x_maximum];

    half_int = maximum(y_smo)/2
    idx_1 = findmin(abs.(y_first_portion .- half_int))
    idx_2 = findmin(abs.(y_second_portion .- half_int))

    hwhm = (x_2[idx_2[2]] .- x_1[idx_1[2]]) ./ 2
    position = x_maximum[1]
    intensity = maximum(y_smo)
    centroid = sum(y_smo ./ sum(y_smo) .* x)

    if y_smo_out == true
        return intensity, position, hwhm, centroid, y_smo
    elseif y_smo_out == false
        return intensity, position, hwhm, centroid
    else
        error("Set y_smo_out to true or false.")
    end
end

"""
    centroid(x::Vector{Float64}, y::Vector{Float64}; smoothing::Bool = false, kwargs...)
    centroid(spectrum::Matrix{Float64}; smoothing::Bool = false, kwargs...)
    centroid(x::Vector{Float64}, y::Matrix{Float64}; smoothing::Bool = false, kwargs...)
    centroid(spectra::Vector{<:Matrix}; smoothing::Bool = false, kwargs...)

Calculate the centroid of a spectrum or a set of spectra.

The function supports optional smoothing and handles various input formats, including single spectra, matrices, and lists of spectra.

# Methods
1. `centroid(x::Vector{Float64}, y::Vector{Float64}; smoothing::Bool = false, kwargs...)`
   Computes the centroid for a single spectrum.
2. `centroid(spectrum::Matrix{Float64}; smoothing::Bool = false, kwargs...)`
   Computes the centroid for a two-column matrix where the first column represents x values and the second column represents y values.
3. `centroid(x::Vector{Float64}, y::Matrix{Float64}; smoothing::Bool = false, kwargs...)`
   Computes the centroid for multiple spectra stored as columns in a matrix.
4. `centroid(spectra::Vector{<:Matrix}; smoothing::Bool = false, kwargs...)`
   Computes the centroids for a list of x-y spectra.

# Arguments
- `x::Vector{Float64}`: The x-axis values (e.g., wavelengths or time points).
- `y::Union{Vector{Float64}, Matrix{Float64}}`: The y-axis values (e.g., intensities), either as a single vector or multiple columns (for multiple spectra).
- `spectrum::Matrix{Float64}`: A two-column matrix with x values in the first column and y values in the second column.
- `spectra::Vector{<:Matrix}`: A list of matrices, each representing an x-y spectrum.
- `smoothing::Bool`: Whether to apply smoothing to the y-axis values before calculating the centroid. Default is `false`.
- `kwargs...`: Additional keyword arguments passed to the smoothing function if smoothing is enabled.

# Returns
- For single spectrum (`x`, `y`): A scalar value representing the centroid position.
- For multiple spectra (`x`, `y` as matrix): A vector where each element is the centroid of a corresponding column in `y`.
- For two-column spectrum (`spectrum`): A scalar value representing the centroid position.
- For list of spectra (`spectra`): A vector where each element is the centroid of a corresponding spectrum in the list.

# Examples
## Example 1: Centroid of a single spectrum
```julia
x = collect(0.:1.:100.)
y_peak = gaussian(x, 1., 50., 10.)
centroid_peak = centroid(x, y_peak)
```

## Example 2: Centroid of a single spectrum with smoothing
```julia
x = collect(0.:1.:100.)
y_peak = gaussian(x, 1., 50., 10.)
centroid_peak = centroid(x, y_peak, smoothing=true, method="gcvspline")
```

## Example 3: Centroid of a single spectrum as an array
```julia
my_spectrum = [x y_peak]
centroid_peak = centroid(my_spectrum)
```

## Example 4: Centroids of an array of y spectra
```julia
ys = [y_peak y_peak y_peak y_peak]
centroid_peaks = centroid(x, ys)
```

## Example 5: Centroids of a vector of x-y spectra
```julia
vector_spectra = [[x y_peak], [x y_peak], [x y_peak]]
centroid_peaks = centroid(vector_spectra)
```

# Notes
- The optional smoothing operation uses an external function (`smooth`) and accepts additional parameters via `kwargs...`.
- Ensure that input dimensions are consistent (e.g., matching lengths for `x` and `y`, or proper formatting for matrices).
"""
function centroid(x::Vector{Float64}, y::Vector{Float64}; smoothing::Bool=false, kwargs...)
    smoothing == true ? out = smooth(x, y; kwargs...) : out = deepcopy(y)
    return sum(out ./ sum(out) .* x)
end
function centroid(spectrum::Matrix{Float64}; smoothing::Bool=false, kwargs...)
    # treat a 2 column spectrum
    return centroid(spectrum[:, 1], spectrum[:, 2]; smoothing=smoothing, kwargs...)
end
function centroid(x::Vector{Float64}, y::Matrix{Float64}; smoothing::Bool=false, kwargs...)
    # treat an array of y spectra
    return [centroid(x, y[:, i]; smoothing=smoothing, kwargs...) for i in 1:size(y, 2)]
end
function centroid(spectra::Vector{<:Matrix}; smoothing::Bool=false, kwargs...)
    # treat a list of x-y spectra
    return [centroid(sp; smoothing=smoothing, kwargs...) for sp in spectra]
end

"""
    find_peaks(
    x::Vector{Float64}, 
    y::Vector{Float64};
    smoothing::Bool = false,
    window_size::Int = 20,
    min_height::Float64 = 0., max_height::Float64 = Inf,
    min_prom::Float64 = 0., max_prom::Float64 = Inf,
    min_width::Float64 = 0., max_width::Float64 = Inf, relwidth::Float64 = 0.5,
    kwargs...
    )

Identify peaks in spectral data (`x`, `y`) based on specified criteria such as height, prominence, and width. 

It optionally applies smoothing to the signal before peak detection.

# Arguments
- `x::Vector{Float64}`: The x-axis values (e.g., wavelengths or time points).
- `y::Vector{Float64}`: The y-axis values (e.g., intensities).
- `smoothing::Bool`: Whether to apply smoothing to the signal before peak detection. Default is `false`.
- `window_size::Int`: Size of the window for peak detection (in indices). Default is `20`.
- `min_height::Float64`: Minimum height of peaks to consider. Default is `0.0`.
- `max_height::Float64`: Maximum height of peaks to consider. Default is `100.0`.
- `min_prom::Float64`: Minimum prominence of peaks to consider. Default is `0.0`.
- `max_prom::Float64`: Maximum prominence of peaks to consider. Default is infinity (`Inf`).
- `min_width::Float64`: Minimum width of peaks to consider (in indices). Default is `0.0`.
- `max_width::Float64`: Maximum width of peaks to consider (in indices). Default is infinity (`Inf`).
- `relwidth::Float64`: Relative width parameter for peak widths. Default is `0.5`.
- Additional keyword arguments (`kwargs...`) are passed to the smoothing function if smoothing is enabled.

# Returns
A named tuple containing:
- `peak_positions::Vector{Float64}`: X positions of detected peaks.
- `peak_heights::Vector{Float64}`: Heights of detected peaks.
- `peak_prominences::Vector{Float64}`: Prominences of detected peaks.
- `peak_hwhms::Vector{Float64}`: Half-width at half maximum values for detected peaks.
- `peak_centroids::Vector{Float64}`: Centroid positions for detected peaks.
- A plot object showing the detected peaks (`plot_peaks`), using the recipe from the Peaks.jl package.

# Examples

## Example 1: Basic peak detection
```julia
x = collect(1:1.0:100)
y = gaussian(x, 1.0, 30.0, 3.0) + lorentzian(x, 1.0, 60.0, 6.0) + 0.01*randn(length(x))

result = find_peaks(x, y; smoothing=false, window_size=10, min_height=0.2)

println("Peak positions: ", result.peak_positions)
println("Peak heights: ", result.peak_heights)
display(result.plot_peaks)
```

## Example 2: Peak detection with gcvspline smoothing
```julia
x = collect(1:1.0:100)
y = gaussian(x, 1.0, 30.0, 3.0) + lorentzian(x, 1.0, 60.0, 6.0) + 0.01*randn(length(x))

result = find_peaks(x, y; smoothing=true, method="gcvspline", window_size=10, min_height=0.2)

println("Peak positions: ", result.peak_positions)
println("Peak heights: ", result.peak_heights)
display(result.plot_peaks)
```

# Notes
- this function is a convenient wrapper around the Peaks.jl functionalities to find peaks in a vector y. For more control, you can directly use the Peaks.jl functions, the package should be available in your environment as it is a dependency for Spectra!
- for smoothing, try the `method = "gcvspline"` or `method = "savgol"` options, which are the most efficient for peak detection.
- the `window_size` parameter is not in the units of `x`, but in the units of index. For example, if you have 1000 points in your x-axis and you want to detect peaks with a window of 10 points, set `window_size=10`.
- the `min_width` and `max_width` parameters are also not in the units of `x`, but in the units of index. For example, if you have 1000 points in your x-axis and you want to detect peaks with a width of 10 points, set `min_width=10` and `max_width=10`.
- the `relwidth` parameter is a relative width parameter for peak widths. It is used to calculate the width of the peak at a certain relative height. For example, if you set `relwidth=0.5`, it will calculate the width of the peak at half its maximum height.

"""
function find_peaks(
    x::Vector{Float64},
    y::Vector{Float64};
    smoothing::Bool=false, # do you want smoothing?
    window_size::Int=20, # what is the size of the window for peak detection? warning => not in the units of X but in the units of index!
    min_height::Float64=0.0,
    max_height::Float64=Inf, # options for peak heights
    min_prom::Float64=0.0,
    max_prom::Float64=Inf, # options for peak prominences
    min_width::Float64=0.0,
    max_width::Float64=Inf,
    relwidth::Float64=0.5, # option for widths, warning min/max are not in the units of X!
    kwargs..., # additional arguments for the smoothing function
)
    # Validate inputs
    if length(x) != length(y)
        error("Vectors 'x' and 'y' must have the same length.")
    end
    if window_size <= 0
        error("'window_size' must be a positive integer.")
    end

    # Apply optional smoothing
    y_input = smoothing ? smooth(x, y; kwargs...) : deepcopy(y)

    # Peak detection steps
    pks_init = findmaxima(y_input, window_size)
    pks_1 = peakheights(pks_init; min=min_height, max=max_height)
    pks_2 = peakproms(pks_1; min=min_prom, max=max_prom)
    pks_3 = peakwidths(pks_2; min=min_width, max=max_width)
    pks_3bis = peakwidths(pks_2; relheight=relwidth)

    # Calculate hwhms
    hwhms = [(x[round(Int, i[2])] - x[round(Int, i[1])]) / 2 for i in pks_3.edges]

    # Calculate centroids
    centroids = Float64[]
    for i in pks_3bis.edges
        lb = x[round(Int, i[1])]
        hb = x[round(Int, i[2])]
        push!(centroids, centroid(x[lb .< x .< hb], y_input[lb .< x .< hb]))
    end

    # Return results as a named tuple
    return (
        peak_positions=x[pks_3.indices],
        peak_heights=pks_3.heights,
        peak_prominences=pks_3.proms,
        peak_hwhms=hwhms,
        peak_centroids=centroids,
        plot_peaks=plotpeaks(
            x, y_input; peaks=pks_3.indices, prominences=true, widths=true
        ),
    )
end

"""
    area_peaks(peak_type::Symbol, amplitude::Float64, hwhm::Float64; lorentzian_fraction=nothing, exponent=nothing)

Calculate the area under Gaussian, Lorentzian, Pseudo-Voigt, or Pearson VII peaks based on their parameters. Areas are calculated using analytical formulas.

# Arguments
- `peak_type::Symbol`: The type of peak. Supported types are:
    - `:gaussian`: Gaussian peak.
    - `:lorentzian`: Lorentzian peak.
    - `:pseudovoigt`: Pseudo-Voigt peak (weighted combination of Gaussian and Lorentzian).
    - `:pearson7`: Pearson VII peak.
- `amplitude::Float64`: Amplitude of the peak (maximum height).
- `hwhm::Float64`: Half-width at half-maximum of the peak.
- `lorentzian_fraction::Union{Float64, Nothing}`: Lorentzian fraction (for Pseudo-Voigt peaks). Must be in [0, 1]. Default is `nothing`.
- `exponent::Union{Float64, Nothing}`: Shape parameter for Pearson VII peaks. Default is `nothing`.

# Returns
- `area::Float64`: The calculated area under the specified peak.

# Examples
## Gaussian Peak
```julia
A = 2.0
hwhm = 0.5
area_gaussian = peak_area("gaussian", amplitude=A, hwhm=hwhm)
```
## Lorentzian Peak
```julia
A = 2.0
hwhm = 0.5
area_lorentzian = peak_area("lorentzian", amplitude=A, hwhm=hwhm)
```
## Pseudo-Voigt Peak
```julia
A = 2.0
hwhm = 0.5
lorentzian_fraction = 0.5
area_pseudovoigt = peak_area("pseudovoigt", amplitude=A, hwhm=hwhm, lorentzian_fraction=lorentzian_fraction)
```
## Pearson VII Peak
```julia
A = 2.0
hwhm = 0.5
exponent = 2.0
area_pearson7 = peak_area("pearson7", amplitude=A, hwhm=hwhm, exponent=exponent)
```
# Notes
1. **Gaussian Formula**:
    `` \\text{Area} = A \\cdot \\text{hwhm} \\cdot \\sqrt{\\frac{\\pi}{\\ln 2}} ``
2. **Lorentzian Formula**:
    `` \\text{Area} = \\pi \\cdot A \\cdot \\text{hwhm} ``
3. **Pseudo-Voigt Formula**:
    `` \\text{Area} = \\eta \\cdot (\\pi \\cdot A \\cdot \\text{hwhm}) + (1-\\eta) \\cdot (A \\cdot \\text{hwhm} \\cdot \\sqrt{\\frac{\\pi}{\\ln 2}}) ``
4. **Pearson VII Formula**:
    `` \\text{Area} = A \\cdot \\text{hwhm} \\cdot \\sqrt{\\frac{\\pi}{2^{1/exponent} - 1}} \\cdot \\frac{\\Gamma(exponent - 0.5)}{\\Gamma(exponent)} ``

# Errors
- Throws an error if an unsupported `peak_type` is provided.
- Throws an error if required parameters for a specific peak type are not provided.
- Throws an error if lorentzian_fraction is not comprised in the [0, 1] interval
"""
function area_peaks(
    peak_type::Symbol,
    amplitude::Float64,
    hwhm::Float64;
    lorentzian_fraction=nothing,
    exponent=nothing,
)
    if peak_type == :gaussian
        return amplitude * hwhm * sqrt(π / log(2))
    elseif peak_type == :lorentzian
        return π * amplitude * hwhm
    elseif peak_type == :pseudovoigt
        lorentzian_fraction = if isnothing(lorentzian_fraction)
            error("lorentzian_fraction required")
        else
            lorentzian_fraction
        end
        if lorentzian_fraction > 1 || lorentzian_fraction < 0
            error("lorentzian_fraction must be in [0, 1]")
        end
        area_L = π * amplitude * hwhm
        area_G = amplitude * hwhm * sqrt(π / log(2))
        return lorentzian_fraction * area_L + (1 - lorentzian_fraction) * area_G
    elseif peak_type == :pearson7
        exponent = isnothing(exponent) ? error("exponent required") : exponent
        k = 2^(1/exponent) - 1
        return amplitude * hwhm * sqrt(π / k) * gamma(exponent - 0.5) / gamma(exponent)
    else
        error("Unsupported peak type")
    end
end
function area_peaks(
    peak_type::Symbol,
    amplitude::Measurement{<:Real},
    hwhm::Measurement{<:Real};
    lorentzian_fraction=nothing,
    exponent=nothing,
)
    if peak_type == :gaussian
        return amplitude * hwhm * sqrt(π / log(2))
    elseif peak_type == :lorentzian
        return π * amplitude * hwhm
    elseif peak_type == :pseudovoigt
        lorentzian_fraction = if isnothing(lorentzian_fraction)
            error("lorentzian_fraction required")
        else
            lorentzian_fraction
        end
        if lorentzian_fraction > 1 || lorentzian_fraction < 0
            error("lorentzian_fraction must be in [0, 1]")
        end
        area_L = π * amplitude * hwhm
        area_G = amplitude * hwhm * sqrt(π / log(2))
        return lorentzian_fraction * area_L + (1 - lorentzian_fraction) * area_G
    elseif peak_type == :pearson7
        exponent = isnothing(exponent) ? error("exponent required") : exponent
        k = 2^(1/exponent) - 1
        return amplitude * hwhm * sqrt(π / k) * gamma(exponent - 0.5) / gamma(exponent)
    else
        error("Unsupported peak type")
    end
end