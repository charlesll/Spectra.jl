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
# baseline.jl contains several "baseline" functions. It is directly dependent on functions.jl.
#
#############################################################################
"""

    baseline(x_input::Vector{Float64}, y_input::Vector{Float64}; roi::Union{Matrix{Float64}, Nothing} = nothing, method::String = "polynomial", kwargs...)
    baseline(x_input::Vector{Float64}, y_input::Matrix{Float64}; roi::Union{Matrix{Float64}, Nothing} = nothing,  method::String = "polynomial", kwargs...)
    baseline(multiple_spectra::Vector{<:Matrix{Float64}}; roi::Union{Matrix{Float64}, Nothing} = nothing,  method::String = "polynomial", kwargs...)

Subtract a baseline from a spectrum `y` sampled at `x` values, or a set of spectra, using specified methods 
that can either fit regions of interests (roi) defined by the user, or use automatic algorithms.

# Arguments
- `x_input::Vector{Float64}`: The x-axis values (e.g., wavelengths or time points).
- `y_input::Union{Vector{Float64}, Matrix{Float64}}`: The spectral data to correct. Can be a single spectrum (vector) or multiple spectra (matrix).
- `multiple_spectra::Vector{<:Matrix{Float64}}`: A collection of spectra, where each spectrum is a matrix with two columns:
    - First column: x-coordinates
    - Second column: y-values (signal intensities)
- `roi::Union{Matrix{Float64}, Nothing}`: Regions of interest for baseline fitting, specified as a matrix where each row defines a range `[start, end]`. Default is `nothing`.
- `method::String`: The baseline fitting method. Default is `"polynomial"`. Supported methods:
    - `"polynomial"` or `"poly"`: Polynomial fitting.
    - `"Dspline"`: 1D spline fitting using Dierckx library.
    - `"gcvspline"`: Generalized Cross Validation Spline fitting.
    - `"exp"`: Exponential background fitting.
    - `"log"`: Logarithmic background fitting.
    - `"rubberband"`: Rubberband baseline fitting.
    - `"als"`: Asymmetric least squares baseline fitting.
    - `"arPLS"`: Asymmetrically reweighted penalized least squares fitting.
    - `"drPLS"`: Doubly reweighted penalized least squares fitting.

# Keyword Arguments
Optional parameters for specific methods:
- `polynomial_order::Int`: Degree of polynomial for polynomial fitting. Default is `1`.
- `s::Float64`: Smoothing coefficient for spline methods. Default is `2.0`.
- `lambda::Float64`: Smoothness parameter for ALS, arPLS, drPLS, and Whittaker methods.
- `p::Float64`: Parameter for ALS algorithm. Default is `0.01`.
- `ratio::Float64`: Parameter for arPLS and drPLS algorithms. Default is `0.01`.
- `niter::Int`: Number of iterations for ALS and drPLS algorithms. Default is `10`.
- `eta::Float64`: Roughness parameter for drPLS algorithm. Default is `0.5`.
- `d::Int`: Order of differences for smoothing algorithms. Default is `2`.
- `d_gcv::Int`: Order of differences for GCV spline algorithm. Default is `3`.
- `p0_exp::Vector{Float64}`: Initial parameters for exponential fitting. Default is `[1., 1., 1.]`.
- `p0_log::Vector{Float64}`: Initial parameters for logarithmic fitting. Default is `[1., 1., 1., 1.]`.

# Returns
- `corrected_signal::Union{Vector{Float64}, Matrix{Float64}}`: Baseline-corrected signal(s).
- `baseline::Union{Vector{Float64}, Matrix{Float64}}`: Calculated baseline(s).

!!! note

    If providing a collection of spectra, you will receive a collection of tuples (corrected_signal, baseline)

# Examples
## Correcting a single spectrum:

```julia
x = collect(50:1.0:500)
background = 10.0 .* sin.(x./50.0) + 0.1.*x
y = 50.0 .* exp.(-log(2) .* ((x .-250.0)./1.0).^2) + background
y_corrected, y_baseline = baseline(x, y, method="drPLS")
```

## Correcting with regions of interest, GCV spline method:
```julia
roi = [[50.0, 200.0], [300.0, 500.0]]
y_corrected, y_baseline = baseline(x, y, roi=roi, method="gcvspline")
```

You can also adjust manually the smoothing spline coefficient `s`:
```julia
y_corrected, y_baseline = baseline(x, y, roi=roi, method="gcvspline", s=1.0)
```

## Using a vector or arrays of y values:
```julia
using Spectra, Plots

# we create a fake signal with 2 peaks plus 2 backgrounds
x = collect(0.:0.2:10.)

# create the signals with create_peaks()
peak_infos = [
    Dict(:type => :gaussian, :amplitude => 10.0, :center => 4.0, :hwhm => 0.6),
    Dict(:type => :gaussian, :amplitude => 5.0, :center => 6.0, :hwhm => 0.4),
]
ys, _ = create_peaks(x, peak_infos)

# add backgrounds
ys[:,1] = ys[:,1] .+ 0.1 * x
ys[:,2] = ys[:,2] .+ 0.2 * x

# fit the background on the first peak: provide as a vector
y_corr, base_ = baseline(x, ys[:,1], method="als")
p1 = plot(x, ys[:,1])
plot!(x, base_)
display(p1)

# Fit the background on multiple peaks: just provide the array!
y_corr, base_ = baseline(x, ys, method="als")
p2 = plot(x, ys, label=["signal 1" "signal 2"])
plot!(x, base_, label=["background 1" "background 2"])
display(p2)
````
## Treating a collection of y values:

Given x1, y1 and x2, y2 two signals, we can create a collection of spectra
and pass it to baseline():
```julia
collection_sp = [[x1 y1], [x2 y2]]
collected_baselines = baseline(collection_sp, method="als)
````
The `collected_baselines` is then a vector that contains tuples of (y_corrected, baseline).
"""
function baseline(
    x_input::Vector{Float64},
    y_input::Vector{Float64};
    roi::Union{Matrix{Float64},Nothing}=nothing,  # Default to `nothing`
    method::String="polynomial",
    kwargs...,
    )

    # Signal standardization
    X_scaler = StatsBase.fit(ZScoreTransform, x_input)
    Y_scaler = StatsBase.fit(ZScoreTransform, y_input)

    # transform the original data
    x = StatsBase.transform(X_scaler, x_input)
    y = StatsBase.transform(Y_scaler, y_input)

    # List of methods that require roi
    methods_requiring_roi = [
        "polynomial", "poly", "Dspline", "gcvspline", "exp", "log", "whittaker"
    ]

    if method in methods_requiring_roi
        # Raise an error if roi is not provided and the method requires it
        if roi === nothing
            error("The 'roi' argument is required for the '$method' method.")
        end

        # Get signals in the roi
        interest_x_unscaled, interest_y_unscaled, interest_index = extract_signal(
            x_input, y_input, roi
        )

        # transform the roi data
        interest_x = StatsBase.transform(X_scaler, interest_x_unscaled)
        interest_y = StatsBase.transform(Y_scaler, interest_y_unscaled)
    end

    if method == "polynomial" || method == "poly"
        # Optional parameters
        poly_order = get(kwargs, :polynomial_order, 1)

        # Fit polynomial
        p = Polynomials.fit(interest_x, interest_y, poly_order)
        baseline_fitted = p.(x)

    elseif method == "Dspline" || method == "unispline"
        # Optional parameters
        splinesmooth = get(kwargs, :s, 2.0)

        # Fit spline using Dierckx
        spl = Spline1D(interest_x, interest_y; k=3, bc="extrapolate", s=splinesmooth)
        baseline_fitted = evaluate(spl, x)

    elseif method == "gcvspline"
        # Optional parameters
        splinesmooth = get(kwargs, :s, nothing)
        d_gcv = get(kwargs, :d_gcv, 3)

        # Use GCVSPL.jl for spline fitting
        # weights = sqrt.(abs.(interest_y))
        # c, wk, ier = Gcvspl.gcv_spl(interest_x, interest_y, weights, splinesmooth)

        # # Evaluate spline at x points
        # baseline_fitted = Gcvspl.spl_der(x, interest_x, c)
        if splinesmooth == nothing
            A = RegularizationSmooth(
                interest_y,
                interest_x,
                d_gcv;
                alg=:gcv_svd,
                extrapolation=ExtrapolationType.Extension,
            )
        else
            A = RegularizationSmooth(
                interest_y,
                interest_x,
                d_gcv;
                λ=splinesmooth,
                alg=:fixed,
                extrapolation=ExtrapolationType.Extension,
            )
        end
        baseline_fitted = A.(x)

    elseif method == "gaussian"
        # Optional parameters
        p0_gauss = get(kwargs, :p0_gaussian, [1.0, 1.0, 1.0])

        # Fit Gaussian
        result = optimize(
            p -> sum((gaussian(interest_x, p...) .- interest_y) .^ 2), p0_gauss
        )
        coeffs = Optim.minimizer(result)

        baseline_fitted = gaussian(x, coeffs...)

    elseif method == "exp"
        # Optional parameters
        p0_exp = get(kwargs, :p0_exp, [1.0, 1.0, 0.0])

        # Fit exponential
        result = optimize(p -> sum((funexp(interest_x, p...) .- interest_y) .^ 2), p0_exp)
        coeffs = Optim.minimizer(result)

        baseline_fitted = funexp(x, coeffs...)

    elseif method == "log"
        # Optional parameters
        p0_log = get(kwargs, :p0_log, [1.0, 1.0, 1.0, 1.0])

        # Fit logarithmic
        result = optimize(p -> sum((funlog(interest_x, p...) .- interest_y) .^ 2), p0_log)
        coeffs = Optim.minimizer(result)

        baseline_fitted = funlog(x, coeffs...)

    elseif method == "rubberband"
        baseline_fitted = rubberband_baseline(x, y)

    elseif method == "whittaker"
        lambda = get(kwargs, :lambda, 1.0e5)

        w = zeros(size(y, 1))
        w[interest_index] .= 1.0

        baseline_fitted = whittaker(x, y, w, lambda)

    elseif method == "als"
        # call ALS algorithm with optional parameters 
        baseline_fitted = als_baseline(
            x,
            y;
            lambda=get(kwargs, :lambda, 1.0e5),
            p=get(kwargs, :p, 0.01),
            niter=get(kwargs, :niter, 10),
        )

    elseif method == "arPLS"
        # call arPLS algorithm with optional parameters
        baseline_fitted = arPLS_baseline(
            x,
            y;
            lambda=get(kwargs, :lambda, 1.0e5),
            ratio=get(kwargs, :ratio, 0.01),
            d=get(kwargs, :d, 2),
        )

    elseif method == "drPLS"
        baseline_fitted = drPLS_baseline(
            x,
            y;
            lambda=get(kwargs, :lambda, 1.0e5),
            niter=get(kwargs, :niter, 10),
            eta=get(kwargs, :eta, 0.5),
            ratio=get(kwargs, :ratio, 0.01),
            d=get(kwargs, :d, 2),
        )

    else
        error("Method not found, check you entered the right name.")
    end

    # Transform back to original scale
    baseline_transformed = vec(StatsBase.reconstruct(Y_scaler, baseline_fitted))
    corrected_signal = y_input .- baseline_transformed

    return corrected_signal, baseline_transformed
end
function baseline(
    x_input::Vector{Float64},
    y_input::Matrix{Float64};
    roi::Union{Matrix{Float64},Nothing}=nothing,  # Default to `nothing`
    method::String="polynomial",
    kwargs...,
    )

    # Initialize an empty matrix for the smoothed output
    y_out = zeros(eltype(y_input), size(y_input))
    base_out = zeros(eltype(y_input), size(y_input))

    # Loop through each column of y and apply the smoothing function
    for i in 1:size(y_input, 2)
        y_out[:, i], base_out[:, i] = baseline(
            x_input, vec(y_input[:, i]); roi=roi, method=method, kwargs...
        )
    end
    return y_out, base_out
end
function baseline(multiple_spectra::Vector{<:Matrix{Float64}};
    roi::Union{Matrix{Float64},Nothing}=nothing,  # Default to `nothing`
    method::String="polynomial",
    kwargs...,
    )
    return [baseline(spectrum[:, 1], spectrum[:, 2]; roi=roi, method=method) for spectrum in multiple_spectra]
end

"""
    als_baseline(x::Vector{Float64}, y::Vector{Float64}; lambda::Float64=1.0e5, p::Float64=0.01, niter::Int=50, d::Int=2) -> Vector{Float64}

Estimate the baseline of a signal using the Asymmetric Least Squares (ALS) method.

# Arguments
- `x::Vector{Float64}`: The x-axis values (must be increasing).
- `y::Vector{Float64}`: The corresponding y-axis values.
- `lambda::Float64=1.0e5`: Smoothing parameter; larger values result in smoother baselines.
- `p::Float64=0.01`: Asymmetry parameter; typically between 0.001 and 0.1.
- `niter::Int=50`: Number of iterations for the algorithm.
- `d::Int=2`: Order of differences for the penalty term.

# Returns
- `z::Vector{Float64}`: The estimated baseline.

!!! note

    This method uses an iterative approach to minimize the asymmetric least squares error, making it suitable for signals with varying background intensity.

# References
G. Eilers and H. Boelens, "Baseline correction with asymmetric least squares smoothing" 2005.
"""
function als_baseline(
    x::Vector{Float64},
    y::Vector{Float64};
    lambda::Float64=1.0e5,
    p::Float64=0.01,
    niter::Int=50,
    d::Int=2,
)

    # Check if x values are equally spaced
    differences = diff(x)
    is_equally_spaced = all(abs.(differences .- differences[1]) .< 1e-10)  # Tolerance for floating-point comparison

    m = length(y)
    E = sparse(1.0I, m, m)  # Identity matrix as a sparse matrix
    w = ones(m)  # Initialize weights to 1

    # Construct the difference matrix
    if is_equally_spaced
        # Compute higher-order differences manually for equally spaced x
        D = E
        for _ in 1:d
            D = diff(D; dims=1)  # Apply diff iteratively to compute higher-order differences
        end
    else
        # Use ddmat for arbitrary spacing
        D = ddmat(x, d)
    end

    # Initialize z with y (baseline estimate)
    z = copy(y)
    for it in 1:niter
        # Update weights matrix W
        W = spdiagm(0 => w)

        # Solve the linear system for updated baseline z
        z = (W + lambda * D'D) \ (w .* y)

        # update weights
        w .= p * (y .> z) .+ (1 - p) * (y .< z)
    end

    return z # Return the estimated baseline
end

"""
    arPLS_baseline(x::Vector{Float64}, y::Vector{Float64}; lambda::Float64=1.0e5, ratio::Float64=0.01, d::Int=2) -> Vector{Float64}

Estimate the baseline of a signal using the Adaptive Reweighted Penalized Least Squares (arPLS) method.

# Arguments
- `x::Vector{Float64}`: The x-axis values (must be increasing).
- `y::Vector{Float64}`: The corresponding y-axis values.
- `lambda::Float64=1.0e5`: Smoothing parameter; larger values result in smoother baselines.
- `ratio::Float64=0.01`: Convergence ratio for stopping criterion.
- `d::Int=2`: Order of differences for the penalty term.

# Returns
- `z::Vector{Float64}`: The estimated baseline.

# Notes
The arPLS algorithm iteratively adjusts weights based on residuals to improve baseline estimation accuracy.

# References
Baek et al., "Baseline correction using adaptive iteratively reweighted penalized least squares," Analyst 140 (2015): 250–257.
"""
function arPLS_baseline(
    x::Vector{Float64}, y::Vector{Float64}; lambda::Float64=1.0e5, ratio::Float64=0.01, d=2
)

    # Check if x values are equally spaced
    differences = diff(x)
    is_equally_spaced = all(abs.(differences .- differences[1]) .< 1e-10)  # Tolerance for floating-point comparison

    m = length(y)
    E = sparse(1.0I, m, m)  # Identity matrix as a sparse matrix
    w = ones(m)  # Initialize weights to 1

    # Construct the difference matrix
    if is_equally_spaced
        # Compute higher-order differences manually for equally spaced x
        D = E
        for _ in 1:d
            D = diff(D; dims=1)  # Apply diff iteratively to compute higher-order differences
        end
    else
        # Use ddmat for arbitrary spacing
        D = ddmat(x, d)
    end

    z = copy(y)  # Initialize z with y (baseline estimate)

    while true
        # Update weights matrix W
        W = spdiagm(0 => w)

        # Solve the linear system for updated baseline z
        z = (W + lambda * D'D) \ (w .* y)

        # Calculate residuals and update weights
        d_residuals = y - z
        dn = d_residuals[d_residuals .< 0]
        m_dn = mean(dn)
        s_dn = std(dn)
        wt = 1.0 ./ (1.0 .+ exp.(2 .* (d_residuals .- (2.0 .* s_dn .- m_dn)) ./ s_dn))

        # check exit condition and backup
        if norm(w-wt)/norm(w) .< ratio
            break
        end

        # Update weights for the next iteration
        w = wt
    end

    return z  # Return the estimated baseline
end

"""
    drPLS_baseline(x::Vector{Float64}, y::Vector{Float64}; lambda::Float64=1.0e5, niter::Int=50, eta::Float64=0.5, ratio::Float64=0.001, d::Int=2) -> Vector{Float64}

Perform baseline correction using the Doubly Reweighted Penalized Least Squares (drPLS) algorithm.

# Arguments
- `x::Vector{Float64}`: The x-axis values (must be increasing).
- `y::Vector{Float64}`: The corresponding y-axis values.
- `lambda::Float64=1.0e5`: Smoothing parameter; larger values result in smoother baselines.
- `niter::Int=50`: Maximum number of iterations for the algorithm.
- `eta::Float64=0.5`: Parameter controlling the influence of weights in the penalty term.
- `ratio::Float64=0.001`: Convergence threshold for stopping criterion.
- `d::Int=2`: Order of differences for the penalty term.

# Returns
- `baseline_fitted::Vector{Float64}`: The estimated baseline.

# Notes
The drPLS algorithm builds on arPLS by introducing an additional reweighting scheme to improve robustness against noise and outliers.

# References
Xu et al., "Baseline correction method based on doubly reweighted penalized least squares," Applied Optics 58 (2019): 3913–3920.
"""
function drPLS_baseline(
    x::Vector{Float64},
    y::Vector{Float64};
    lambda::Float64=1.0e5,
    niter::Int=50,
    eta::Float64=0.5,
    ratio::Float64=0.01,
    d::Int=2,
)
    # Check if x values are equally spaced
    differences = diff(x)
    is_equally_spaced = all(abs.(differences .- differences[1]) .< 1e-10)  # Tolerance for floating-point comparison

    m = length(y)
    E = sparse(1.0I, m, m)

    # Create sparse identity matrix
    I_n = sparse(1.0I, m, m)

    # Create sparse difference matrices
    if is_equally_spaced
        # Compute higher-order differences manually for equally spaced x
        D = E
        for _ in 1:d
            D = diff(D; dims=1)  # Apply diff iteratively to compute higher-order differences
        end
        D_1 = diff(I_n; dims=1)
    else
        # Use ddmat for arbitrary spacing
        D = ddmat(x, d)
        D_1 = ddmat(x, 1)
    end

    # final construction
    D = D'D  # Penalty matrix for smoothness
    D_1 = D_1'D_1  # First-order difference matrix

    # Initialize weights and other variables
    w_0 = ones(m)
    w = copy(w_0)
    Z = copy(w_0)

    # Fitting procedure
    for jj in 1:niter
        W = spdiagm(0 => w)  # Diagonal weight matrix
        Z_prev = copy(Z)

        # Update baseline estimate
        Z = (W + D_1 + lambda * (I_n - eta * W) * D) \ (W * y)

        # Check convergence
        if norm(Z - Z_prev) / norm(Z_prev) < ratio
            break
        end

        # Update weights based on residuals
        d_residuals = y - Z
        d_negative = d_residuals[d_residuals .< 0]
        sigma_negative = std(d_negative)
        mean_negative = mean(d_negative)

        w .=
            0.5 * (
                1 .-
                exp.(
                    jj * (d_residuals .- (-mean_negative + 2 * sigma_negative)) ./
                    sigma_negative,
                ) ./ (
                    1 .+ abs.(
                        exp.(
                            jj * (d_residuals .- (-mean_negative + 2 * sigma_negative)) ./
                            sigma_negative,
                        ),
                    )
                )
            )
    end

    return Z  # Return the fitted baseline
end

"""
    rubberband_baseline(x::AbstractVector, y::AbstractVector) -> baseline::Vector

Compute the baseline of a spectrum using the rubberband (convex hull) method.

# Arguments
- `x::AbstractVector` : 1D array of the independent variable (e.g., wavenumber or wavelength). Can be unequally spaced but must be sorted or sortable.
- `y::AbstractVector` : 1D array of the dependent variable (e.g., absorbance or intensity), same length as `x`.

# Returns
- `baseline::Vector` : The estimated baseline of `y`, obtained by constructing a convex hull underneath the spectrum and linearly interpolating between hull points.

# Description
This function estimates the baseline of a spectrum by stretching a "rubberband" under the data points. It first identifies the lower convex hull of the `(x, y)` points, then performs linear interpolation between consecutive hull points to produce the baseline. Subtracting this baseline from the original spectrum can be used for baseline-corrected spectra.

# Example
```julia
x = 1:100
y = 0.05 .* x .+ sin.(0.1 .* x) .+ 0.5 .* exp.(-0.01 .* (x .- 50).^2)
baseline = rubberband_baseline(x, y)
y_corrected = y .- baseline
```
"""
function rubberband_baseline(x::AbstractVector, y::AbstractVector)
    if length(x) != length(y)
        error("x and y must have the same length.")
    end
    # Ensure x and y are sorted by x
    if !issorted(x)
        idx = sortperm(x)
        x = x[idx]
        y = y[idx]
    end

    n = length(x)
    baseline = zeros(n)

    # Start convex hull from the first point
    hull_indices = [1]

    for i in 2:n
        push!(hull_indices, i)
        while length(hull_indices) >= 3
            j = length(hull_indices)
            # Three points: (x1,y1),(x2,y2),(x3,y3)
            x1, y1 = x[hull_indices[j-2]], y[hull_indices[j-2]]
            x2, y2 = x[hull_indices[j-1]], y[hull_indices[j-1]]
            x3, y3 = x[hull_indices[j]], y[hull_indices[j]]
            # Compute cross product to check convexity
            if (x2 - x1)*(y3 - y1) - (y2 - y1)*(x3 - x1) <= 0
                deleteat!(hull_indices, j-1)  # remove middle point
            else
                break
            end
        end
    end

    # Interpolate baseline along hull points
    for i in 1:length(hull_indices)-1
        i1, i2 = hull_indices[i], hull_indices[i+1]
        baseline[i1:i2] .= y[i1] .+ (y[i2]-y[i1]) .* ((x[i1:i2] .- x[i1]) ./ (x[i2]-x[i1]))
    end

    return baseline
end
