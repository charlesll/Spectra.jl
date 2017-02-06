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

function poly(p::Vector{Float64},x::Array{Float64})
    segments = zeros(size(x)[1],size(p)[1])
    for i = 1:size(p)[1]
        segments[:,i] = p[i].*x[:].^(i-1)
    end
    return sum(segments,2)
end

# From https://rosettacode.org/wiki/Polynomial_regression#Julia
function polyfit(x::Array{Float64}, y::Array{Float64}, n::Int64)
  A = [ float(x[i])^p for i = 1:length(x), p = 0:n ]
  return A \ y
end

"""
For Gaussian peaks in spectra
	gaussiennes(amplitude::Array{Float64},centre::Array{Float64},hwhm::Array{Float64},x::Array{Float64},style::String = "None")
"""
function gaussiennes(amplitude::Array{Float64},centre::Array{Float64},hwhm::Array{Float64},x::Array{Float64};style::String = "None")
    segments = zeros(size(x)[1],size(amplitude)[1])
    if style == "None"
        for i = 1:size(amplitude)[1]
            segments[:,i] = amplitude[i] .*exp(-log(2) .* ((x[:,1]-centre[i])./hwhm[i]).^2)
        end
    elseif style == "poly"
        for i = 1:size(amplitude)[1]
            segments[:,i] = poly(vec(amplitude[i,:]),x[:,2]) .*exp(-log(2) .* ((x[:,1]-(poly(vec(centre[i,:]),x[:,2])))./poly(vec(hwhm[i,:]),x[:,2])).^2)
        end	    	
    else
        error("Not implemented, see documentation")
    end
    return sum(segments,2), segments
end

function lorentziennes(amplitude::Array{Float64},centre::Array{Float64},hwhm::Array{Float64},x::Array{Float64};style::String = "None")
    segments = zeros(size(x)[1],size(amplitude)[1])
    if style == "None"
        for i = 1:size(amplitude)[1]
            segments[:,i] = amplitude[i] ./ (1 + ((x[:,1]-centre[i])./hwhm[i]).^2)
        end
    elseif style == "poly"
        for i = 1:size(amplitude)[1]
            segments[:,i] = poly(vec(amplitude[i,:]),x[:,2]) ./ (1 + ((x[:,1]-(poly(vec(centre[i,:]),x[:,2])))./poly(vec(hwhm[i,:]),x[:,2])).^2)
        end	    	
    else
        error("Not implemented, see documentation")
    end
    return sum(segments,2), segments
end

function pearson7(a1::Array{Float64},a2::Array{Float64},a3::Array{Float64},a4::Array{Float64},x::Array{Float64};style::String = "None")
    segments = zeros(size(x)[1],size(a1)[1])
    if style == "None"
        for i = 1:size(a1)[1]
            segments[:,i] = a1[i] ./ (1 + ((x[:,1]-a2[i])./a3[i]).^2 .* (2.0.^(1./a4[i]) - 1.0))
        end
    elseif style == "poly"
        for i = 1:size(a1)[1]
            segments[:,i] = poly(vec(a1[i,:]),x[:,2]) ./ (1 + ((x[:,1]-(poly(vec(a2[i,:]),x[:,2])))./poly(vec(a3[i,:]),x[:,2])).^2 .* (2.0.^(1./vec(a4[i,:])) - 1.0))
        end	    	
    else
        error("Not implemented, see documentation")
    end
    return sum(segments,2), segments
end

function pseudovoigts(amplitude::Array{Float64},centre::Array{Float64},hwhm::Array{Float64},lorentzian_fraction::Array{Float64},x::Array{Float64};style::String = "None")
    segments = zeros(size(x)[1],size(amplitude)[1])
	test1 = find(lorentzian_fraction .<0.0)
	test2 = find(lorentzian_fraction .>1.0)
    if size(test1)[1] != 0 || size(test2)[1] != 0
		error("lorentzian_fraction should be comprised between 0 and 1")
	end
	if style == "None"
		for i = 1:size(amplitude)[1]
			segments[:,i] =  amplitude[i].*( (1.0 - lorentzian_fraction[i]) .*exp(-log(2) .* ((x[:,1]-centre[i])./hwhm[i]).^2)   + lorentzian_fraction[i] ./ (1 + ((x[:,1]-centre[i])./hwhm[i]).^2))
        end
    elseif style == "poly"
        for i = 1:size(amplitude)[1]
			segments[:,i] =  vec(lorentzian_fraction[i,:]) .* (poly(vec(amplitude[i,:]),x[:,2]) ./ (1 + ((x[:,1]-(poly(vec(centre[i,:]),x[:,2])))./poly(vec(hwhm[i,:]),x[:,2])).^2))     +      (1- vec(lorentzian_fraction[i,:])) .* poly(vec(amplitude[i,:]),x[:,2]) .*exp(-log(2) .* ((x[:,1]-(poly(vec(centre[i,:]),x[:,2])))./poly(vec(hwhm[i,:]),x[:,2])).^2) 
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

# to correct for shifts in X between two spectra, for the inversion
function xshift_inversion(data::Array{Float64},p::Array{Float64})
    xaxis = data[:,1]
    shifted1 = data[:,2]

    spl = Spline1D(xaxis-p[1], shifted1.*p[2] + shifted1.^2.*p[3])
    y = evaluate(spl,xaxis)
end

# to correct a spectrum for a p shift in X
function xshift_direct(original_x::Array{Float64}, original_y::Array{Float64}, p::Float64)
    spl = Spline1D(original_x-p, original_y)
    corrected_y = evaluate(spl,original_x)
	return original_x, corrected_y, p
end

# to correct a shift between two spectra using a reference peak
function xshift_correction(full_x::Array{Float64}, full_shifted_y::Array{Float64}, ref_x::Array{Float64}, ref_y::Array{Float64},shifted_y::Array{Float64})
	fit = curve_fit(xshift_inversion, [ref_x shifted_y], ref_y, [1.0, 1.0,1.0])
	parameters = fit.param
	return xshift_direct(full_x, full_shifted_y, parameters[1])
end

"""
Savitzky-Golay filter of window half-width M and degree N
M is the number of points before and after to interpolate, i.e. the full width
of the window is 2M+1
Code from https://medium.com/@acidflask/smoothing-data-with-julia-s-generated-functions-c80e240e05f3#.45v03x6it
"""
immutable SavitzkyGolayFilter{M,N} end

@generated function Base.call{M,N,T}(::Type{SavitzkyGolayFilter{M,N}},
     data::AbstractVector{T})

     #Create Jacobian matrix
     J = zeros(2M+1, N+1)
     for i=1:2M+1, j=1:N+1
         J[i, j] = (i-M-1)^(j-1)
     end
     e₁ = zeros(N+1) 
     e₁[1] = 1.0
    
     #Compute filter coefficients
     C = J' \ e₁

     #Evaluate filter on data matrix

     To = typeof(C[1] * one(T)) #Calculate type of output
     expr = quote
         n = size(data, 1)
         smoothed = zeros($To, n)
         @inbounds for i in eachindex(smoothed)
             smoothed[i] += $(C[M+1])*data[i]
         end
         smoothed
     end

     for j=1:M
         insert!(expr.args[6].args[2].args[2].args, 1,
             :(if i - $j ≥ 1
                 smoothed[i] += $(C[M+1-j])*data[i-$j]
               end)
         )
         push!(expr.args[6].args[2].args[2].args,
             :(if i + $j ≤ n
                 smoothed[i] += $(C[M+1+j])*data[i+$j]
               end)
         )
     end

     return expr
end