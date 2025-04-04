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
# Helper functions for baseline methods
function gaussian(x, a, b, c)
    return a * exp.(-log(2) * ((x .- b) ./ c).^2)
end

function funexp(x, a, b, c)
    return a * exp.(b * (x .- c))
end

function funlog(x, a, b, c, d)
    return a * log.(-b * (x .- c)) .- d * x.^2
end

"""
	poly(p::Vector{Float64},x::Vector{Float64})

Builds a polynomial curve given parameters `p::Vector{Float64}` at the `x::Vector{Float64}` values.

For a linear curve, p = [1.0,1.0], for a second order polynomial, p = [1.0,1.0,1.0], etc.;

"""
function poly(p::Vector{Float64}, x::Vector{Float64})
    segments = zeros(size(x)[1],size(p)[1])
    for i = 1:size(p)[1]
        segments[:,i] = p[i].*x.^(i-1)
    end
    return sum(segments, dims=2)
end

"""
	gaussiennes(amplitude::Array{Float64},centre::Array{Float64},hwhm::Array{Float64},x::Array{Float64};style::String = "None")

gaussiennes, written in the plural french form there, is a function that allows to build gaussian peaks. The gaussian function used there is:

    y = amplitude x exp(-ln(2) x [(x-centre)/hwhm]^2 )

You can enter the amplitude, centre and half-width at half-maximum (hwhm) values as arrays of float 64 (even containing one float value), without specifying style. hwhm is proportional to the standard deviation sigma:

	hwhm= sqrt(2xln(2)) x sigma

that is used in a normal distribution (see function normal_dist).

# Inputs

	amplitude: Array{Float64}
		peaks amplitudes
	centre: Array{Float64}
		peaks centres
	hwhm: Array{Float64}
		peaks half-width at middle heights (hwhm);
	x: Array{Float64}
		x axis values;

# Options

	style: ASCIIString = "None"
		see examples below.

# Outputs

	y_calc: Array{Float64}
		calculated y values
	y_peaks: Array{Float64}
		calculated y values of the different peaks.

# Examples

To have four gaussian peaks centered at 800, 900, 1000 and 1100 cm-1 with hwhm of 50 cm-1 on a Raman spectrum, you will enter:

	```julia-repl
	julia> y_calc, y_peaks = gaussiennes([1.0,1.0,1.0,1.0], [800.0,900.0,1000.0,1100.0], [50.0,50.0,50.0,50.0], x)
    ```

and y_peaks will contain in 4 columns the 4 different y values of the peaks, and y_calc their sum (the total model). Now, if you want to calculate more complex models, such as for instance contructing how the Raman peaks of water vary with pressure, you might like to parametrize the variations of the peak parameters rather than just fitting each spectrum. This will provide more robust fits of the spectra, as you will fit them together, and will also force you to find the correct underlying mathematical assumption.

The gaussiennes function allows you to do that. If you specify style = "poly", you can enter arrays for the amplitudes, centres and half-widths at half-maximum (hwhm) of the peaks, with in each column the coefficients for the polynomial variations of this parameters. The second column of x will need to contain the second variable for those polynomial functions.

Let's say for instance that we have one peak at 900 cm-1 in a pure material. It's frequency seems to linearly shift with increasing the amount of hydrogen in this material, but it's intensity is non-linearly increasing, following a quadratic variation. It's width seems constant.

How to write that with gaussiennes? Well, first you need to construct a relevant x axis: first column contains the frequency, and the second one contains the chemical variable value. In our case, we want to model the peak between 800 and 1000 cm-1, for 1 wt% H. So we have an x array build like:

	```julia-repl
    julia> frequency = collect(800:1:1000)
    julia> x = ones(length(frequency),2)
    julia> x[:,1] = frequency[:]
    julia> x[:,2] = 1.0
	```

Ok, now lets build our y peaks:

	```julia-repl
    julia> amplitudes = [1.0 0.1 0.1]
    julia> frequencies = [900.0 2.0]
    julia> hwhm = 20.0
    julia> y_calc, y_peaks = gaussiennes(amplitudes, frequencies, hwhm, x)
	```

This should provide you how the shape of the peak is as a function of both the frequency and the chemical composition there. If you want to go further, you might just want to stick gaussiennes in a loop, and play with creating various peaks with changing the chemical parameter in the x[:,2] column!
"""
function gaussiennes(amplitude::Array{Float64},centre::Array{Float64},hwhm::Array{Float64},x::Array{Float64};style::String = "None")
    segments = zeros(size(x)[1],size(amplitude)[1])
    if style == "None"
        for i = 1:size(amplitude)[1]
            segments[:,i] = amplitude[i] .*exp.(-log.(2.0) .* ((x[:,1].-centre[i])./hwhm[i]).^2)
        end
    elseif style == "poly"
        for i = 1:size(amplitude)[1]
            segments[:,i] = poly(vec(amplitude[i,:]),x[:,2]) .*exp.(-log.(2) .* ((x[:,1].-(poly(vec(centre[i,:]),x[:,2])))./poly(vec(hwhm[i,:]),x[:,2])).^2)
        end
    else
        error("Not implemented, see documentation")
    end
    return sum(segments, dims=2), segments
end

"""
	lorentziennes(amplitude::Array{Float64},centre::Array{Float64},hwhm::Array{Float64},x::Array{Float64};style::String = "None")

Inputs
------

	amplitude: Array{Float64}
		peaks amplitudes
	centre: Array{Float64}
		peaks centres
	hwhm: Array{Float64}
		peaks half-width at middle heights (hwhm)
	x: Array{Float64}
		x axis values

Options
-------

	style: ASCIIString = "None", see examples in the gaussiennes documentation.

Outputs
-------

	y_calc: Array{Float64}
		calculated y values
	y_peaks: Array{Float64}
		calculated y values of the different peaks.

"""
function lorentziennes(amplitude::Array{Float64},centre::Array{Float64},hwhm::Array{Float64},x::Array{Float64};style::String = "None")
    segments = zeros(size(x)[1],size(amplitude)[1])
    if style == "None"
        for i = 1:size(amplitude)[1]
            segments[:,i] = amplitude[i] ./ (1.0 .+ ((x[:,1].-centre[i])./hwhm[i]).^2)
        end
    elseif style == "poly"
        for i = 1:size(amplitude)[1]
            segments[:,i] = poly(vec(amplitude[i,:]),x[:,2]) ./ (1.0 .+ ((x[:,1].-(poly(vec(centre[i,:]),x[:,2])))./poly(vec(hwhm[i,:]),x[:,2])).^2)
        end
    else
        error("Not implemented, see documentation")
    end
    return sum(segments, dims=2), segments
end

"""
	pearson7(a1::Array{Float64},a2::Array{Float64},a3::Array{Float64},a4::Array{Float64},x::Array{Float64};style::String = "None")

a Pearson 7 peak with formula a1 ./ (1 + ((x-a2)./a3).^2 .* (2.0.^(1./a4) - 1.0))

Inputs
------

	a1: Array{Float64}
		parameter a1
	a2: Array{Float64}
		parameter a2
	a3: Array{Float64}
		parameter a3
	a4: Array{Float64}
		parameter a4
	x: Array{Float64}
		x axis values

Options
-------

	style: ASCIIString = "None", see examples in the gaussiennes documentation.

Outputs
-------

	y_calc: Array{Float64}
		calculated y values
	y_peaks: Array{Float64}
		y values of the different peaks.

"""
function pearson7(a1::Array{Float64},a2::Array{Float64},a3::Array{Float64},a4::Array{Float64},x::Array{Float64};style::String = "None")
    segments = zeros(size(x)[1],size(a1)[1])
    if style == "None"
        for i = 1:size(a1)[1]
            segments[:,i] = a1[i] ./ (1.0 .+ ((x[:,1].-a2[i]) ./a3[i]) .^2 .* (2.0 .^(1.0 ./a4[i]) .- 1.0))
        end
    elseif style == "poly"
        for i = 1:size(a1)[1]
            segments[:,i] = poly(vec(a1[i,:]),x[:,2]) ./ (1. .+ ((x[:,1].-(poly(vec(a2[i,:]),x[:,2]))) ./poly(vec(a3[i,:]),x[:,2])) .^2 .* (2.0 .^(1,0 ./vec(a4[i,:])) .- 1.0))
        end
    else
        error("Not implemented, see documentation")
    end
    return sum(segments, dims=2), segments
end

"""
	pseudovoigts(amplitude::Array{Float64},centre::Array{Float64},hwhm::Array{Float64},lorentzian_fraction::Array{Float64},x::Array{Float64};style::String = "None")

A mixture of gaussian and lorentzian peaks.

Inputs
------
	amplitude: Array{Float64}
		peaks amplitudes
	centre: Array{Float64}
		peaks centres
	hwhm: Array{Float64}
		peaks half-width at middle heights (hwhm)
	lorentzian_fraction: Array{Float64}
		lorentzian fraction of the pseudovoigt function. Should be comprised between 0 and 1;
	x: Array{Float64}
		x axis values
Options
-------

	style: ASCIIString = "None", see examples in the gaussiennes documentation.

Outputs
-------
	y_calc: Array{Float64}
		calculated y values
	y_peaks: Array{Float64}
		y values of the different peaks

"""
function pseudovoigts(amplitude::Array{Float64},centre::Array{Float64},hwhm::Array{Float64},lorentzian_fraction::Array{Float64},x::Array{Float64};style::String = "None")
    segments = zeros(size(x)[1],size(amplitude)[1])
	test1 = findall(lorentzian_fraction .< 0.0)
	test2 = findall(lorentzian_fraction .> 1.0)
    if size(test1)[1] != 0 || size(test2)[1] != 0
		error("lorentzian_fraction should be comprised between 0 and 1")
	end
	if style == "None"
		for i = 1:size(amplitude)[1]
			segments[:,i] =  amplitude[i].*( (1.0 .- lorentzian_fraction[i]) .*exp.(-log(2) .* ((x[:,1].-centre[i])./hwhm[i]).^2)   .+ lorentzian_fraction[i] ./ (1.0 .+ ((x[:,1].-centre[i])./hwhm[i]).^2))
        end
    elseif style == "poly"
        for i = 1:size(amplitude)[1]
			segments[:,i] =  vec(lorentzian_fraction[i,:]) .* (poly(vec(amplitude[i,:]),x[:,2]) ./ (1 .+ ((x[:,1].-(poly(vec(centre[i,:]),x[:,2])))./poly(vec(hwhm[i,:]),x[:,2])).^2))     .+      (1.0 .- vec(lorentzian_fraction[i,:])) .* poly(vec(amplitude[i,:]),x[:,2]) .*exp.(-log(2) .* ((x[:,1].-(poly(vec(centre[i,:]),x[:,2])))./poly(vec(hwhm[i,:]),x[:,2])).^2)
        end
    else
        error("Not implemented, see documentation")
	end
    return sum(segments, dims=2), segments
end


"""
	normal_dist(nd_amplitudes::Array{Float64},nd_centres::Array{Float64},nd_sigmas::Array{Float64},x::Array{Float64})

The real normal distribution / gaussian function

Inputs
------

	nd_amplitudes: Array{Float64}
		distribution probability
	nd_centres: Array{Float64}
		distribution position
	nd_sigmas: Array{Float64}
		sigma value of the distribution
	x: Array{Float64}
		x axis values

Options
-------

	style: ASCIIString = "None", see examples in the gaussiennes documentation.

Outputs
-------

	y_calc: Array{Float64}
		calculated y values
	y_peaks: Array{Float64}
		y values of the different peaks

"""
function normal_dist(nd_amplitudes::Array{Float64},nd_centres::Array{Float64},nd_sigmas::Array{Float64},x::Array{Float64})
    segments = zeros(size(x)[1],size(nd_amplitudes)[1])
    for i = 1:size(nd_amplitudes)[1]
        segments[:,i] = nd_amplitudes[i]./(nd_sigmas[i].*sqrt(2.0.*pi))  .*  exp.(- (x[:,1].-nd_centres[i]).^2 ./ (2.0.*nd_sigmas[i].^2.0))
    end
    return sum(segments, dims=2), segments
end

