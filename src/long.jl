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
# Long Correction
# Charles Le Losq
# RSES-ANU 2016
# Long's correction of Raman spectra and normalisation
# last rev. Oct 2016 for convertion in Julia
# ensures strictly increasing values of wavenumber
# calc. e.s.e. as Long cor. norm. sqrt(n_raw) 3d output col.
# exp. format to avoid null ese.
# program long3;
# Original program in Pascal written by J. Roux, modified for Julia by C. Le Losq.
#############################################################################

function long(data::Array{Float64},temp::Float64,wave::Float64)

    h::Float64 = 6.62606896e-34   # J.s    Plank constant
    k::Float64 = 1.38066e-23      # J/K    Boltzman
    c::Float64 = 2.9979e8         # m/s    Speed of light
    nu0::Float64 = 1.0./wave*1e9     # nu0 laser is in m-1
    T::Float64 = temp + 273.15    # K temperature

    x::Float64 = data[]:,1]
    y::Float64 = data[:,2]

    # Calculate the error on data as sqrt(y). If y <= 0, then error = abs(y).
    ese::Float64 = sqrt(abs(y))./abs(y) # relative errors

    # then we proceed to the correction (Neuville and Mysen, 1996; Le Losq et al., 2012)
    nu::Float64 = 100.0.*x # cm-1 -> m-1 Raman shift
    t0::Float64 = nu0.^3.*nu./(nu0-nu)
    t1::Float64 = 1 - exp(-h.*c.*nu./(k.*T)) # c in m/s  : t1 dimensionless
    long::Float64 = y.*t0.*t1;% pour les y
    
    long = long./trapz(x,long) # area normalisation

    eselong::Float64 = ese.*long # error calculation
    
    return x, long, eselong

