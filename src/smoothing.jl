
"""
    whittaker(x::Vector{Float64}, y::Vector{Float64}, w::Vector{Float64}, lambda::Float64; d::Int=2) -> Vector{Float64}

Smooth a signal using the Whittaker smoother with divided differences.

# Arguments
- `x::Vector{Float64}`: The x-axis values of the data (must be increasing).
- `y::Vector{Float64}`: The corresponding y-axis values, assumed to be sampled at equal intervals.
- `w::Vector{Float64}`: Weights for each data point. Higher weights indicate greater importance in smoothing.
- `lambda::Float64`: Smoothing parameter; larger values result in smoother outputs. It is recommended to be set based on the length of `x`.
- `d::Int=2`: Order of differences for the penalty term (default is 2).

# Returns
- `z::Vector{Float64}`: The smoothed y-axis values.

# Notes
- For equally spaced x values, higher-order differences are computed manually.
- For unequally spaced x values, a difference matrix is constructed using the `ddmat` function.
- A sparse matrix representation is used for computational efficiency.

# References
- Eilers, P.H.C. (2003). "A perfect smoother." Analytical Chemistry, 75, 3631–3636.

# Example

```julia
using Spectra, Random, Plots
x = sort(rand(50) .* 10)                 # Randomly spaced x values in [0, 10]
true_y = sin.(x)
y = true_y .+ 0.1 .* randn(length(x))  # Noisy sine wave
w = ones(length(x))                     # Equal weights
lambda = 0.1                          # Smoothing parameter, !!! its value depends on the length of x !!!
z = whittaker(x, y, w, lambda; d=2)

p1 = plot(x, true_y); plot!(x, y); plot!(x, z), display(p1)
```

Matlab version by Paul Eilers, 2003
Julia translation by Charles Le Losq 2017, revised 2025
"""
function whittaker(
    x::Vector{Float64}, y::Vector{Float64}, w::Vector{Float64}, lambda::Float64; d::Int=2
)
    if length(x) != length(y)
        error("x must be the same length as y")
    elseif length(y) != length(w)
        error("w must be the same length as y")
    end

    # Check if x values are equally spaced
    differences = diff(x)
    is_equally_spaced = all(abs.(differences .- differences[1]) .< 1e-10)  # Tolerance for floating-point comparison

    # Smoothing
    m = length(y)
    E = sparse(1.0I, m, m)  # Identity matrix as a sparse matrix

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

    # Construct the diagonal matrix for weights
    W = spdiagm(0 => w)  # Diagonal sparse matrix

    # Solving the linear system
    z = (W + lambda * D'D) \ (w .* y)

    return z
end

"""

    ddmat(x::Vector{Float64}, d::Int)

Construct the difference matrix for arbitrary spacing of a vector x, with d the order of differences

# Arguments
- `x::Vector{Float64}`: sampling positions
- `d::Int`: order of differences

# Output
- `D::SparseMatrixCSC{Float64, Int}`: the difference matrix

# Notes
- `D * Y`: gives divided differences of order d
- Matlab version in Paul Eilers, 2003, Julia translation by Charles Le Losq 2017, revised 2025
"""
function ddmat(x::Vector{Float64}, d::Int)
    m = length(x)
    if d == 0
        return sparse(1.0I, m, m) # Sparse identity matrix
    else
        dx = x[(d + 1):m] .- x[1:(m - d)]
        V = spdiagm(m - d, m - d, 0 => 1.0 ./ dx)
        D_prev = ddmat(x, d - 1)     # Recursive call
        D_diff = diff(D_prev; dims=1) # Apply diff along rows
        return V * D_diff[1:(m - d), :] # Ensure dimensions match for multiplication
    end
end

"""
    smooth(x::Vector{Float64}, y::Union{Vector{Float64}, Matrix{Float64}}; 
           method::String = "gcvspline", 
           window_length::Int = 5, 
           polyorder::Int = 2, 
           lambda::Float64 = 10.0^5, 
           d::Int = 2, 
           w = nothing, 
           d_gcv::Int = 3)

Smooth a signal using various smoothing methods.

# Arguments
- `x::Vector{Float64}`: The x-coordinates of the input signal (e.g., time or spatial values).
- `y::Union{Vector{Float64}, Matrix{Float64}}`: The y-coordinates of the input signal to be smoothed.
- `method::String`: The smoothing method to use. Options include:
    - `"whittaker"`: Whittaker smoother, which uses weights (`w`) and smoothing parameter (`lambda`).
    - `"gcvspline"`: Generalized cross-validation spline smoothing (requires `DataInterpolations.jl`).
    - `"flat"`: Moving average.
    - `"hanning"`, `"hamming"`, `"bartlett"`, `"blackman"`: Window-based smoothing methods (requires `DSP.jl`).
    - `"savgol"`: Savitzky-Golay filter.
- `window_length::Int`: The length of the smoothing window (used in window-based and Savitzky-Golay methods). Must be a positive odd integer
- `polyorder::Int`: The polynomial order for the Savitzky-Golay filter. Must be less than `window_length`
- `lambda::Float64`: Smoothing parameter for the Whittaker smoother. Higher values result in smoother fits
- `d::Int`: Order of differences for the Whittaker smoother. Default is 2.
- `w`: Weights for the Whittaker smoother. size(w) should be equal to size(y). If not provided (`nothing`), uniform weights are used by default.
- `d_gcv::Int`: Order of differences for the GCV spline algorithm (used in `"gcvspline"`)

# Returns
- A smoothed vector or matrix of signals (`Vector{Float64}` or `Matrix{Float64}`).

# Notes
- **Deprecated Methods**: The methods `"GCVSmoothedNSpline"`, `"MSESmoothedNSpline"`, and `"DOFSmoothedNSpline"` are no longer supported and will throw an error if used.
- For window-based methods (`"flat"`, `"hanning"`, etc.), the signal is symmetrically extended at both ends to reduce edge effects.
- The Savitzky-Golay filter requires that `window_length` be a positive odd integer and that `polyorder` be less than `window_length`.
- For Whittaker smoothing, if weights (`w`) are provided, they must have the same length as `x`.

# Example

```julia
using Spectra, Plots

x = collect(1:10)
y = [1.0, 2.5, 3.0, 4.2, 5.1, 6.0, 7.3, 8.1, 9.4, 10.0]
method = "hanning"
window_length = 3

smoothed_y = smooth(x, y; method=method, window_length=window_length)
p1 = plot(x, [y smoothed_y], label=["Original" "Smoothed"], title="Smoothing Example", xlabel="X-axis", ylabel="Y-axis")
display(p1)
```

# Errors
- Throws an error if an unsupported or deprecated method is specified.
- Throws an error if `window_length` is not a positive odd integer.
- For Savitzky-Golay smoothing:
    - Throws an error if `polyorder` is greater than or equal to `window_length`.
- For Whittaker smoothing:
    - Throws an error if weights (`w`) are provided but do not match the length of `x`.

# References
- For Window-based filtering,see DSP.jl documentation: [https://docs.juliadsp.org/stable/windows/](https://docs.juliadsp.org/stable/windows/)
"""
function smooth(
    x::Vector{Float64},
    y::Vector{Float64};
    method::String="gcvspline",
    window_length::Int=5,
    polyorder::Int=2,
    lambda::Float64=10.0^5,
    d::Int=2,
    w=nothing,
    d_gcv::Int=3,
)
    isodd(window_length) || error("window_length must be a positive odd integer")
    window_length > 0 || error("window_length must be a positive odd integer")
    polyorder > 0 || error("polyorder must be a positive integer")
    polyorder < window_length || error("polyorder must be less than window_length")

    if method in ["GCVSmoothedNSpline", "MSESmoothedNSpline", "DOFSmoothedNSpline"]
        error(
            "The methods GCVSmoothedNSpline, MSESmoothedNSpline, and DOFSmoothedNSpline are deprecated and no longer supported.",
        )

    elseif method == "gcvspline"
        # from the DataInterpolations package
        smoother = RegularizationSmooth(y, x, d_gcv; alg=:gcv_svd, extrapolation = ExtrapolationType.Constant)
        return smoother.(x)

    elseif method == "whittaker"
        if w == nothing
            w = ones(size(y))
        end
        return whittaker(x, y, w, lambda; d=d)

    elseif method in ["flat", "hanning", "hamming", "bartlett", "blackman"]
        # We use the DSP package, which allows a nearly direct translation of my Python code (thanks Perplexity AI)

        # Step 1: Extend the signal symmetrically
        s = vcat(
            reverse(y[2:window_length]), y, reverse(y[(end - window_length + 1):(end - 1)])
        )

        # Step 2: Create the smoothing window
        w = if method == "flat"  # Moving average
            ones(window_length)
        else
            getfield(DSP, Symbol(method))(window_length)  # Generate window using DSP methods
        end

        # Step 3: Normalize the window and apply convolution
        y_filt = conv(w ./ sum(w), s)

        # Step 4: Adjust for edge effects
        shift = Int((length(y_filt) - length(y)) / 2)
        return y_filt[(shift + 1):(end - shift)]

    elseif method == "savgol"
        # Savitzky-Golay filter
        # OLD CODE WITH STAGEDFILTERS, CREATE ISSUES AT THE BOUNDARY
        # Create type-stable output vector
        #smoothed = zeros(eltype(y), length(y)); # <--- wholesome, type stable code.
        #return smooth!(SavitzkyGolayFilter{window_length,polyorder}, y, smoothed)

        # NEW CODE WITH SavitzkyGolay.jl, no problem at the boundaries
        sg = savitzky_golay(y, window_length, polyorder)
        return sg.y
        
    else
        error("Unknown smoothing method: $method")
    end
end
function smooth(
    x::Vector{Float64},
    y::Matrix{Float64};
    method::String="whittaker",
    window_length::Int=5,
    polyorder::Int=2,
    lambda::Float64=10.0^5,
    d::Int=2,
    w=nothing,
    d_gcv::Int=3,
)
    # Check if y is a matrix
    if ndims(y) != 2
        error("y must be a matrix")
    end

    # initialize weights in advance
    if w == nothing
        w = ones(size(y))
    else # check the size of w compared to y
        if size(w) != size(y)
            error("w must be the same size as y")
        end
    end

    # Initialize an empty matrix for the smoothed output
    smoothed = zeros(eltype(y), size(y))

    # Loop through each column of y and apply the smoothing function
    for i in 1:size(y, 2)
        smoothed[:, i] = smooth(
            x,
            vec(y[:, i]);
            method=method,
            window_length=window_length,
            polyorder=polyorder,
            lambda=lambda,
            d=d,
            w=w,
            d_gcv=d_gcv,
        )
    end
    return smoothed
end
