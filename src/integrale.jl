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
# This files contains functions to perform calculations of numeric integrale values on real dataset x y.
#
#
#############################################################################

"""
The 'trapz' function allows to perform numerical integration of a signal.
	integrale (float64) = trapz(x:vector of float 64, y:vector of float 64)
"""
function trapz{Tx<:Number, Ty<:Number}(x::Vector{Tx}, y::Vector{Ty})
    # Trapezoidal integration rule
    local n = length(x)
    if (length(y) != n)
        error("Vectors 'x', 'y' must be of same length")
    end
    r = zero(zero(Tx) + zero(Ty))
    if n == 1; return r; end
    for i in 2:n
        r += (x[i] - x[i-1]) * (y[i] + y[i-1])
    end
    trapz_int = r/2
    return trapz_int
end

"""
    function gaussianarea(Amplitude::Array{Float64},HWHM::Array{Float64}; eseAmplitude::Array{Float64} = 0, eseHWHM::Array{Float64} = 0)
Return the integrated itnensity (area) of a gaussian with inputed amplitude and HWHM.
Options are eseAmplitude (0 by default) and eseHWHM (0 by default).

If Amplitude and HWHM are arrays, outputs will be arrays.
"""
function gaussianarea(Amplitude::Array{Float64},HWHM::Array{Float64}; eseAmplitude::Array{Float64} = [0.0], eseHWHM::Array{Float64} = [0.0])
    
    area::Array{Float64} = sqrt(pi./log(2)).*Amplitude.*HWHM
    if (eseAmplitude == 0.0) && (eseHWHM == 0.0)
        esearea::Array{Float64} = 0.0
    else
        esearea = sqrt((pi./log(2).*HWHM).^2 .* eseAmplitude.^2 + (pi./log(2).*Amplitude).^2 .* eseHWHM.^2)
    end
    return area, esearea
end