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
# Temperature and laser wavelength effects correction for Raman spectra
# Charles Le Losq
# RSES-ANU 2016
# Long's correction of Raman spectra and normalisation
# last rev. Oct 2016 for convertion in Julia
# ensures strictly increasing values of wavenumber
# calc. e.s.e. as Long cor. norm. sqrt(n_raw) 3d output col.
# exp. format to avoid null ese.
#
# THOSE FORMULA ARE WRITTEN FOR STOKES RAMAN. May be easily adapted for anti-Stokes upon request.
#
# References
# Shuker, Reuben, and Robert Gammon. 1970. “Raman-Scattering Selection-Rule Breaking and the Density of States in Amorphous Materials.” Physical Review Letters 25 (4): 222–25.
# Galeener, F. L., and Sen, P. N. 1978. “Theory of the First-Order Vibrational Spectra of Disordered Solids.” Physical Review B 17 (4): 1928–33.
# Mysen, B. O., L. W. Finger, D. Virgo, and F. A. Seifert. 1982. “Curve-Fitting of Raman Spectra of Silicate Glasses.” American Mineralogist 67: 686–95.
# Brooker et al. 1988 Assessment of correction procedures for reduction of Raman spectra. Journal of Raman Spectroscopy 19(2), 71-78.
# Neuville, D. R., and B. O. Mysen. 1996. “Role of Aluminium in the Silicate Network: In Situ, High-Temperature Study of Glasses and Melts on the Join SiO₂-NaAl0₂.” Geochimica et Cosmochimica Acta 60: 1727–37.
# Le Losq, C., D. R. Neuville, R. Moretti, and J. Roux. 2012. “Determination of Water Content in Silicate Glasses Using Raman Spectrometry: Implications for the Study of Explosive Volcanism.” American Mineralogist 97 (5-6): 779–90. doi:10.2138/am.2012.3831.
# Hehlen, B. 2010. “Inter-Tetrahedra Bond Angle of Permanently Densified Silicas Extracted from Their Raman Spectra.” Journal of Physics: Condensed Matter 22 (2): 025401.
#############################################################################

function tlcorrection(data::Array{Float64},temp::Float64,wave::Float64;correction="long",density=2210.0)

    h::Float64 = 6.626070040e-34   # J S    Plank constant from NIST
	hb::Float64 = 1.054571800e-34 # J S    Reduced Plank constant from NIST
    k::Float64 = 1.38064852e-23      # J K-1    Boltzman constant from NIST
    c::Float64 = 299792458.         # M S-1    Speed of light from NIST
    nu0::Float64 = 1.0./wave*1e9     # nu0 laser is in M-1 (wave is in nm)
    T::Float64 = temp + 273.15    # K temperature
	# density is in KG M-3

    x::Array{Float64} = data[:,1]
    y::Array{Float64} = data[:,2]

    # Calculate the error on data as sqrt(y). If y <= 0, then error = abs(y).
    ese::Array{Float64} = sqrt(abs(y))./abs(y) # relative errors

	# get the Raman shift in m-1
	nu::Array{Float64} = 100.0.*x # cm-1 -> m-1 Raman shift
	
    # then we proceed to the correction 
	if correction == "long"
		# Formula used in Mysen et al. (1982), Neuville and Mysen (1996) and Le Losq et al. (2012) (corrected for using the Planck constant in the last reference)
		# It is that reported in Brooker et al. (1988) with the addition of a scaling nu0^3 coefficient for adimentionality
    	frequency::Array{Float64} = nu0.^3.*nu./((nu0-nu).^4) # frequency correction; dimensionless
    	boltzman::Array{Float64} = 1 - exp(-h.*c.*nu./(k.*T)) # temperature correction with Boltzman distribution; dimensionless
    	ycorr::Array{Float64} = y.*frequency.*boltzman; # correction
		
	elseif correction == "galeener"
		# This uses the formula reported in Galeener and Sen (1978) and Brooker et al. (1988); it uses the Bose-Einstein / Boltzman distribution
		# Formula from  without the scaling vo^3 coefficient reported in Mysen et al. (1982), Neuville and Mysen (1996) and Le Losq et al. (2012)
    	frequency = nu./((nu0-nu).^4) # frequency correction; M^3
    	boltzman = 1 - exp(-h.*c.*nu./(k.*T)) # temperature correction with Boltzman distribution; dimensionless
    	ycorr = y.*frequency.*boltzman; # correction
		
	elseif correction =="hehlen"
		# this uses the formula reported in Hehlen et al. 2010
    	frequency = nu./(nu0.^3.*density) # frequency + density correction; M/KG
    	boltzman = 1 - exp(-h.*c.*nu./(k.*T)) # dimensionless
    	ycorr = y.*frequency.*boltzman; # correction
	else
		error("Not implemented, choose between long, galeener or hehlen.")
	end
	
    ycorr = ycorr./trapz(x,ycorr) # area normalisation

    esecorr::Array{Float64} = ese.*ycorr # error calculation

    return x, ycorr, esecorr
end
