#############################################################################
#Copyright (c) 2016-2019 Charles Le Losq
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

    peakmeas(x::Array{Float64,1}, y::Array{Float64,1}; smoothing = "yes", method = "savgol", window_length=5, polyorder=2, ese_y=1., y_smo_out=false)

The peakmeas function allows performing measurements of the position, width, intensity and centroïd of a dominant peak in a provided x-y signal.

It smooths the signal with a Savitzky-Golay filter prior to measuring the peak position, width and intensity. It is advised to check that the M and N values of the Savitzky-Golay filter are adequate for your problem before trusting the results from peakmeas. For that, just use the y_smo_out option.

half-width at half-maximum are calculated as the width of the peak at half its maximum intensity. This calculation is not affected by any asumption of peak symmetry (no fitting is done).

Centroïd is calculated as sum(y./sum(y).*x).

INPUTS:

	x: Array{Float64}, the x values;

	y: Array{Float64}, the y values.

OPTIONS:

	smoothing, String, triggers the smoothing of the spectrum if set to yes (default value);

	filter, Symbol, the filter that will be used. See the smooth function documentation;

	M=5, the M parameter for smoothing y with a Savitzky-Golay filter. See smooth function documentation;

	N=2, the M parameter for smoothing y with a Savitzky-Golay filter. See smooth function documentation;

	y_smo_out=false, the smoothed signal. Signal will be smoothed if set to true, using the SavitzkyGolayFilter function with the M and N values. y_smo output will also be provided.

OUTPUTS:

	intensity: Float64, the intensity of the peak;

	position: Float64, the position of the peak;

	hwhm: Float64, the half-width at half-maximum of the peak;

	centroïd: Float64, the centroïd of the peak;

	if y_smo_out is set to true, then another output is provided:

	y_smo: Array{Float64}, the smoothed y signal.
"""
function peakmeas(x::Array{Float64,1}, y::Array{Float64,1}; smoothing = "yes", method = "savgol", window_length=5, polyorder=2, ese_y=1., y_smo_out=false)
    ### PRELIMINARY CHECK: INCREASING SIGNAL
    if x[end,1] < x[1,1]
        x = flipdim(x,1)
        y = flipdim(y,1)
    end

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
