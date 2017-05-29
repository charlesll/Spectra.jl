#############################################################################
#Copyright (c) 2016 Charles Le Losq
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
	peakhw(x::Array{Float64},y::Array{Float64};M=5,N=2,y_smo_out=false)

The peakhw function allows performing measurements of the position, width and intensity of a peak.

It also allows smoothing the signal with a Savitzky-Golay filter prior to measuring the peak position, width and intensity, see the options.

It is advised to check that the M and N values of the Savitzky-Golay filter are adequate for your problem before trusting the results from peakhw. For that, just use the y_smo_out option. 

INPUTS:
	
	x: Array{Float64}, the x values;

	y: Array{Float64}, the y values.

OPTIONS:

	M=5, the M parameter for smoothing y with a Savitzky-Golay filter. See SavitzkyGolayFilter documentation;

	N=2, the M parameter for smoothing y with a Savitzky-Golay filter. See SavitzkyGolayFilter documentation;

	y_smo_out=false, the smoothed signal. Signal will be smoothed if set to true, using the SavitzkyGolayFilter function with the M and N values. y_smo output will also be provided.
	
OUTPUTS:

	x_maximum: Float64, the position of the peak;
	
	hwhm: Float64, the half-width at half-maximum of the peak;
	
	if y_smo_out is set to true, then another output is provided:
	
	y_smo: Array{Float64}, the smoothed y signal.
"""
function peakhw(x::Array{Float64},y::Array{Float64};M=5,N=2,y_smo_out=false)
    ### PRELIMINARY CHECK: INCREASING SIGNAL
    if x[end,1] < x[1,1]
        x = flipdim(x,1)
    y = flipdim(y,1)
    end

    y_smo = SavitzkyGolayFilter{M,N}(vec(y))

    x_maximum = x[y_smo .== maximum(y_smo)]
    x_1 = x[x .<x_maximum]
    x_2 = x[x .>=x_maximum]
    y_first_portion = y_smo[x .<x_maximum]
    y_second_portion = y_smo[x .>=x_maximum]
    half_int = maximum(y_smo)/2
    idx_1 = findmin(abs(y_first_portion-half_int))
    idx_2 = findmin(abs(y_second_portion-half_int))
    hwhm = (x_2[idx_2[2]]-x_1[idx_1[2]])./2

    if y_smo_out == true
      return x_maximum[1], hwhm, y_smo
    elseif y_smo_out ==false
      return x_maximum[1], hwhm
    else
      error("Set y_smo_out to true or false.")
    end

end
