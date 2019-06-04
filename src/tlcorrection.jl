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
#############################################################################

"""
	tlcorrection(data::Array{Float64},temp::Float64,wave::Float64;correction="long",normalisation="area",density=2210.0)

Inputs
------
	data: Array{Float64}
		input spectrum with x and y in first and second columns respectively
	temp: Float64
		temperature in °C
	wave: Float64
		wavenumber at which the spectrum was acquirred in nm

Options
-------
	correction: String
		equation used for the correction. Choose between "long", "galeener", or "hehlen". Default = "long".
	normalisation: String
		indicate if you want to normalise your signal or not. Choose between "intensity", "area", or "no". Default = "area".
	density: Float64
		density of the studied material in kg m-3, to be used with the "hehlen" equation. Default = 2210.0 (density of silica).

Outputs
-------

	x: Array{Float64}
		x values
	long: Array{Float64}
		corrected y values
	eselong: Array{Float64}
		errors calculated as sqrt(y) on raw data and propagated after the correction.

Notes
-------

This correction uses the formula reported in Galeener and Sen (1978), Mysen et al. (1982), Brooker et al. (1988) and Hehlen et al. (2010).

The "galeener" equation is the exact one reported in Galeener and Sen (1978), which is a modification from Shuker and Gammon (1970) for accounting of (vo - v)^4 dependence of the Raman intensity. See also Brooker et al. (1988) for further discussion.

The "long" equation is that of Galeener and Sen (1978) corrected by a vo^3 coefficient for removing the cubic meter dimension of the equation of "galeener". This equation has been used in Mysen et al. (1982), Neuville and Mysen (1996) and Le Losq et al. (2012).

The "hehlen" equation is that reported in Hehlen et al. (2010). It actually originates before this publication (Brooker et al. 1988). It uses a different correction that avoid crushing the signal below 500 cm-1. THerefore, it has the advantage of keeping intact the Boson peak signal in glasses.
"""
function tlcorrection(data::Array{Float64},temp::Float64,wave::Float64;correction="long",normalisation="area",density=2210.0)

    h::Float64 = 6.626070040e-34   # J S    Plank constant from NIST
	hb::Float64 = 1.054571800e-34 # J S    Reduced Plank constant from NIST
    k::Float64 = 1.38064852e-23      # J K-1    Boltzman constant from NIST
    c::Float64 = 299792458.0         # M S-1    Speed of light from NIST
    nu0::Float64 = 1.0 ./wave .*1.0e9     # nu0 laser is in M-1 (wave is in nm)
    T::Float64 = temp + 273.15    # K temperature
	# density is in KG M-3

    x::Array{Float64} = data[:,1]
    y::Array{Float64} = data[:,2]

    # Calculate the error on data as sqrt(y). If y <= 0, then error = abs(y).
    ese::Array{Float64} = sqrt.(abs.(y))./abs.(y) # relative errors

	# get the Raman shift in m-1
	nu::Array{Float64} = 100.0.*x # cm-1 -> m-1 Raman shift

    # then we proceed to the correction
	if correction == "long"
		# Formula used in Mysen et al. (1982), Neuville and Mysen (1996) and Le Losq et al. (2012) (corrected for using the Planck constant in the last reference)
		# It is that reported in Brooker et al. (1988) with the addition of a scaling nu0^3 coefficient for adimentionality
    	frequency::Array{Float64} = nu0.^3.0.*nu./((nu0.-nu).^4) # frequency correction; dimensionless
    	boltzman::Array{Float64} = 1.0 .- exp.(-h.*c.*nu./(k.*T)) # temperature correction with Boltzman distribution; dimensionless
    	ycorr::Array{Float64} = y.*frequency.*boltzman; # correction

	elseif correction == "galeener"
		# This uses the formula reported in Galeener and Sen (1978) and Brooker et al. (1988); it uses the Bose-Einstein / Boltzman distribution
		# Formula from  without the scaling vo^3 coefficient reported in Mysen et al. (1982), Neuville and Mysen (1996) and Le Losq et al. (2012)
    	frequency = nu./((nu0.-nu).^4) # frequency correction; M^3
    	boltzman = 1.0 .- exp.(-h.*c.*nu./(k.*T)) # temperature correction with Boltzman distribution; dimensionless
    	ycorr = y.*frequency.*boltzman; # correction

	elseif correction =="hehlen"
		# this uses the formula reported in Hehlen et al. 2010
    	frequency = 1.0./(nu0.^3.0.*density) # frequency + density correction; M/KG
    	boltzman = 1.0 .- exp.(-h.*c.*nu./(k.*T)) # dimensionless
    	ycorr = nu.*y.*frequency.*boltzman; # correction
	else
		error("Not implemented, choose between long, galeener or hehlen.")
	end

	if normalisation == "area"
    	ycorr = ycorr./trapz(x,ycorr) # area normalisation
	elseif normalisation == "intensity"
		ycorr= ycorr./maximum(ycorr) # max. intensity normalisation
	elseif normalisation == "no"
		# nothing will happen
	else
		error("Set the optional normalisation parameter to area, intensity or no.")
	end

    esecorr::Array{Float64} = ese.*ycorr # error calculation

    return x, ycorr, esecorr
end
