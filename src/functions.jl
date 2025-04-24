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
	normal_dist(x, amplitude, mean, sigma)
	normal_dist(x, p)

Normal distribution with parameters amplitude, mean, sigma or a vector p = [amplitude, mean, sigma]

# Notes:
- normal_dist(x, p) is a shorthand for normal_dist(x, p[1], p[2], p[3])
"""
function normal_dist(x, amplitude::Float64, mean::Float64, sigma::Float64)
    return amplitude ./ (sigma .* sqrt(2.0 .* pi)) .*
           exp.(- (x .- mean) .^ 2 ./ (2.0 .* sigma .^ 2.0))
end
function normal_dist(x, p::Vector{Float64})
    return normal_dist(x, p[1], p[2], p[3]) # p = [amplitude, mean, sigma]
end

"""
	gaussian(x, amplitude, center, hwhm)
	gaussian(x, p)

Gaussian function with parameters amplitude, center, hwhm or a vector p = [amplitude, center, hwhm]

# Notes:
- the half-width at half maximum (hwhm) of the gaussian peak is related to the standard deviation sigma by: hwhm = sqrt(2*log(2)) * sigma
- gaussian(x, p) is a shorthand for gaussian(x, p[1], p[2], p[3])
"""
function gaussian(x, amplitude, center, hwhm)
    return @. amplitude * exp.(-log(2) * ((x .- center) ./ hwhm) .^ 2)
end
function gaussian(x, p)
    if length(p) > 3
        error("p should contain 3 values: amplitude, center, hwhm")
    end
    return gaussian(x, p[1], p[2], p[3])# p = [amplitude, center, hwhm]
end

"""
	lorentzian(x, amplitude, center, hwhm)
	lorentzian(x, p)

Lorentzian function with parameters amplitude, center, hwhm or a vector p = [amplitude, center, hwhm]

# Notes:
- hwhm: half-width at half maximum
- lorentzian(x, p) is a shorthand for lorentzian(x, p[1], p[2], p[3])
"""
function lorentzian(x, amplitude, center, hwhm)
    return @. amplitude ./ (1.0 .+ ((x .- center) ./ hwhm) .^ 2)
end
function lorentzian(x, p)
    if length(p) > 3
        error("p should contain 3 values: amplitude, center, hwhm")
    end
    # p = [amplitude, center, hwhm]
    return lorentzian(x, p[1], p[2], p[3])
end

"""
	pseudovoigt(x, amplitude, center, hwhm, lorentzian_fraction)
	pseudovoigt(x, p)

Pseudovoigt function with parameters amplitude, center, hwhm, lorentzian_fraction or a vector p = [amplitude, center, hwhm, lorentzian_fraction]

Calculated as lorentzian_fraction*lorentzian + (1 - lorentzian_fraction)*gaussian

# Notes:
- hwhm: half-width at half maximum
- pseudovoigt(x, p) is a shorthand for lorentzian(x, p[1], p[2], p[3], p[4])
- lorentzian_fraction is a value between 0 and 1 that controls the mixing between Gaussian and Lorentzian
"""
function pseudovoigt(x, amplitude, center, hwhm, lorentzian_fraction)
    g = gaussian(x, amplitude, center, hwhm)
    l = lorentzian(x, amplitude, center, hwhm)
    return @. (1-lorentzian_fraction) * g + lorentzian_fraction * l
end
function pseudovoigt(x, p)
    # p = [amplitude, center, hwhm, lorentzian_fraction]
    return pseudovoigt(x, p[1], p[2], p[3], p[4])
end

"""
	pearson7(x, amplitude, center, hwhm, exponent)
	pearson7(x, p)

Pearson 7 function with parameters amplitude, center, hwhm, exponent or a vector p = [amplitude, center, hwhm, exponent]

Equation is amplitude / (1 + ((x - center)/hwhm)^2 * (2^(1/exponent) - 1))^exponent

# Notes:
- pearson7(x, p) is a shorthand for pearson7(x, p[1], p[2], p[3], p[4])
"""
function pearson7(x, amplitude, center, hwhm, exponent)
    return @. amplitude ./
           (1.0 .+ ((x .- center) ./ hwhm) .^ 2 .* (2.0 .^ (1.0 ./ exponent) .- 1.0)) .^
           exponent
end
function pearson7(x, p)
    return pearson7(x, p[1], p[2], p[3], p[4]) # p = [amplitude, center, hwhm, shape]
end

"""
Generates multiple peaks and their sum from a collection of peak descriptions.

# Arguments
- `x::Vector{Float64}`: X-axis values where peaks will be evaluated.
- `peak_infos::Vector{Dict}`: List of dictionaries describing peaks. Each dictionary should contain:
    - `:type`: Peak type (`:gaussian`, `:lorentzian`, `:pseudovoigt`, `:pearson7`)
    - Required parameters for the specified peak type:
        - Gaussian: `:amplitude`, `:center`, `:hwhm`
        - Lorentzian: `:amplitude`, `:center`, `:hwhm`
        - PseudoVoigt: `:amplitude`, `:center`, `:hwhm`, `:lorentzian_fraction`
        - Pearson7: `:amplitude`, `:center`, `:hwhm`, `:exponent`

# Returns
- `peaks::Matrix{Float64}`: Matrix where each column represents a individual peak
- `total_spectrum::Vector{Float64}`: Sum of all peaks

# Examples

```julia
x = collect(0:0.1:10)
peak_infos = [
Dict(:type => :gaussian, :amplitude => 1.0, :center => 5.0, :hwhm => 0.5),
Dict(:type => :lorentzian, :amplitude => 0.8, :center => 7.0, :hwhm => 1.2),
Dict(:type => :pearson7, :amplitude => 0.8, :center => 3.0, :hwhm => 0.2, :exponent => 1.9)
]
peaks, total = create_peaks(x, peak_infos)
```

"""
function create_peaks(x::Vector{Float64}, peak_infos::Vector{Dict{Symbol,Any}})
    # Initialize output matrix
    peaks = zeros(length(x), length(peak_infos))

    # Define peak function dispatch dictionary
    peak_functions = Dict(
        :gaussian => gaussian,
        :lorentzian => lorentzian,
        :pseudovoigt => pseudovoigt,
        :pearson7 => pearson7,
    )

    # Process each peak definition
    for (i, info) in enumerate(peak_infos)
        peak_type = info[:type]

        # Get required parameters based on peak type
        params = if peak_type == :gaussian || peak_type == :lorentzian
            (info[:amplitude], info[:center], info[:hwhm])
        elseif peak_type == :pseudovoigt
            (info[:amplitude], info[:center], info[:hwhm], info[:lorentzian_fraction])
        elseif peak_type == :pearson7
            (info[:amplitude], info[:center], info[:hwhm], info[:exponent])
        else
            error("Unsupported peak type: $peak_type")
        end

        # Generate peak and store in matrix column
        peaks[:, i] = peak_functions[peak_type](x, params...)
    end

    # Calculate total spectrum
    total_spectrum = sum(peaks; dims=2)

    return peaks, vec(total_spectrum)
end

####################
# HELPER FUNCTIONS #
####################
function funexp(x, a, b, c)
    return a * exp.(b * (x .- c))
end

function funlog(x, a, b, c, d)
    return a * log.(-b * (x .- c)) .- d * x .^ 2
end

"""
	poly(p::Vector{Float64},x::Vector{Float64})

Builds a polynomial curve given parameters `p::Vector{Float64}` at the `x::Vector{Float64}` values.

For a linear curve, p = [1.0,1.0], for a second order polynomial, p = [1.0,1.0,1.0], etc.;

"""
function poly(p::Vector{Float64}, x::Vector{Float64})
    segments = zeros(size(x)[1], size(p)[1])
    for i in 1:size(p)[1]
        segments[:, i] = p[i] .* x .^ (i-1)
    end
    return sum(segments; dims=2)
end

