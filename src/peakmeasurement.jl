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

The peakmeas function allows performing measurements of the position, width, intensity and centroid of a dominant peak in a provided `x`-`y` signal.

It smooths the signal with a Savitzky-Golay filter prior to measuring the peak position, width and intensity. It is advised to check that the M and N values of the Savitzky-Golay filter are adequate for your problem before trusting the results from peakmeas. For that, just use the y_smo_out option.

half-width at half-maximum are calculated as the width of the peak at half its maximum intensity. This calculation is not affected by any asumption of peak symmetry (no fitting is done).

Centroïd is calculated as sum(y./sum(y).*x).

# Inputs

	`x::Array{Float64}`: the x values
	`y::Array{Float64}`: the y values

# Options

    `smoothing::String`: triggers the smoothing of the spectrum if set to yes (default value);
    `filter::Symbol`: the filter that will be used. See the smooth function documentation;
    `M::Int`: M parameter for smoothing y with a Savitzky-Golay filter. See smooth function documentation. Default = 5.
    `N::Int`: N parameter for smoothing y with a Savitzky-Golay filter. See smooth function documentation. Default = 2. 
    `y_smo_out::bool`: Outputs the smoothed signal.

# Outputs

	`intensity::Float64`: peak intensity
	`position::Float64`: peak position
	`hwhm: Float6:: peak` half-width at half-maximum
	`centroïd::Float64`: peak centroid
"""
function peakmeas(x::Array{Float64,1}, y::Array{Float64,1}; smoothing = "yes", method = "savgol", window_length=5, polyorder=2, ese_y=1., y_smo_out=false)
    ### PRELIMINARY CHECK: INCREASING SIGNAL
    sp = flipsp([x y])
    x = deepcopy(sp[:,1])
    y = deepcopy(sp[:,2])

	if smoothing == "yes"
    	y_smo = vec(smooth(x,y,method = method, window_length=window_length, polyorder=polyorder))
	else
		y_smo = collect(y)
	end

    x_maximum = x[y_smo .== maximum(y_smo)]

	if size(x_maximum,1) >= 2
		error("Something went wrong, is there two peaks with identical intensities in the signal? This function treats signal with one main peak...")
	end

    x_1 = x[x .<x_maximum];
	x_2 = x[x .>=x_maximum];
    y_first_portion = y_smo[x .<x_maximum];
	y_second_portion = y_smo[x .>=x_maximum];

    half_int = maximum(y_smo)/2
    idx_1 = findmin(abs.(y_first_portion.-half_int))
    idx_2 = findmin(abs.(y_second_portion.-half_int))

    hwhm = (x_2[idx_2[2]].-x_1[idx_1[2]])./2
    position = x_maximum[1]
    intensity = maximum(y_smo)
    centroid = sum(y_smo./sum(y_smo).*x)

    if y_smo_out == true
      return intensity, position, hwhm, centroid, y_smo
    elseif y_smo_out ==false
      return intensity, position, hwhm, centroid
    else
      error("Set y_smo_out to true or false.")
    end

end

"""

    centroid(x::Array{Float64,2}, y::Array{Float64,2}; smoothing::Bool = false, kwargs...)

Returns the centroids of pairs of x (n by m) - y (n by m) signal(s) of n frequencies, m samples. To smooth the signal(s), set `smoothing::Bool` to `true`. See [`smooth`](@ref) for optional smoothing arguments. 

# Notes

Centroids are calculated as ``sum(y[:,i] / sum(y[:,i]) .* x[:,i])`` `

# Example

```julia
Using Spectra
x = collect(0.:1.:100.)
y_sum, ys = gaussiennes([10., 20.],[50., 60.], [3., 2.], x)
centroid(x, y_sum)

# output

1×1 Matrix{Float64}:
55.71428571428571

```

"""
function centroid(x::Array{Float64}, y::Array{Float64}; smoothing::Bool = false, kwargs...)

    y_ = copy(y)

    if smoothing
        for i in 1:size(x, 2)
            y_[:, i] = smooth(x[:, i], y[:, i]; kwargs...)
        end
    end

    return sum(y_ ./ sum(y_, dims=1) .* x, dims=1)
end