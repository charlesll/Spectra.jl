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

# this function smooths a peak and measures its height and its half-width
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
    hwhm = x_2[idx_2[2]]-x_1[idx_1[2]]

    if y_smo_out == true
      return x_maximum, hwhm, y_smo
    elseif y_smo_out ==false
      return x_maximum, hwhm
    else
      error("Set y_smo_out to true or false.")
    end

end
