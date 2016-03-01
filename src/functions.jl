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
Bootstrap function
    bootstrap(data::Array{Float64}, ese::Array{Float64},nbsample::Int64)

Bootstrap of Raman spectra. We generate new datapoints with the basis of existing data and their standard deviation
"""
function bootstrap(data::Array{Float64}, ese::Array{Float64},nbsample::Int64)
    bootsamples = zeros(size(data)[1],nbsample)
    for i = 1:nbsample
        randn!(bootsamples[:,i]) .* ese[:] + data[:]
    end
    return bootsamples
end

function poly(p::Vector{Float64},x::Array{Float64})
    segments = zeros(size(x)[1],size(p)[1])
    for i = 1:size(p)[1]
        segments[:,i] = p[i].*x[:].^(i-1)
    end
    return sum(segments,2)
end
"""
For Gaussian peaks in spectra
	gaussiennes(g_amplitudes::Array{Float64},g_frequency::Array{Float64},g_hwhm::Array{Float64},x::Array{Float64},style::ASCIIString = "None")
"""
function gaussiennes(g_amplitudes::Array{Float64},g_frequency::Array{Float64},g_hwhm::Array{Float64},x::Array{Float64};style::ASCIIString = "None")
    segments = zeros(size(x)[1],size(g_amplitudes)[1])
    if style == "None"
        for i = 1:size(g_amplitudes)[1]
            segments[:,i] = g_amplitudes[i] .*exp(-log(2) .* ((x[:,1]-g_frequency[i])./g_hwhm[i]).^2)
        end
    elseif style == "poly"
        for i = 1:size(g_amplitudes)[1]
            segments[:,i] = poly(squeeze(g_amplitudes[i,:],1),x[:,2]) .*exp(-log(2) .* ((x[:,1]-(poly(squeeze(g_frequency[i,:],1),x[:,2])))./poly(squeeze(g_hwhm[i,:],1),x[:,2])).^2)
        end	    	
    else
        error("Not implemented, see documentation")
    end
    return sum(segments,2), segments
end


#The real normal distribution / gaussian function
function normal_dist(nd_amplitudes::Array{Float64},nd_centres::Array{Float64},nd_sigmas::Array{Float64},x::Array{Float64})
    segments = zeros(size(x)[1],size(nd_amplitudes)[1])
    for i = 1:size(nd_amplitudes)[1]
        segments[:,i] = nd_amplitudes[i]./(nd_sigmas[i].*sqrt(2*pi))  .*  exp(- (x[:,1]-nd_centres[i]).^2 ./ (2.*nd_sigmas[i].^2))
    end
    return sum(segments,2), segments
end