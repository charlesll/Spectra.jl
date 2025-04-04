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

"""
	tlcorrection(spectrum::Matrix{Float64}, temperature_C::Float64, laser_wavelength::Float64;
                 correction="long", normalisation="area", density=2210.0)
    tlcorrection(x::Vector{Float64}, y::Vector{Float64}, temperature_C::Float64, laser_wavelength::Float64;
                 correction="long", normalisation="area", density=2210.0)
    tlcorrection(multiple_spectra::Vector{<:Matrix{Float64}}, temperature_C::Float64, laser_wavelength::Float64;
                 correction="long", normalisation="area", density=2210.0)

Temperature and laser wavelength correction for Raman spectra using one of three 
available correction equations: "long", "galeener", or "hehlen". 
Also supports optional normalization of the corrected spectra.

# Arguments
- `spectrum::Matrix{Float64}`: A single spectrum with x (Raman shift in cm⁻¹) and y (intensity) values in the first and second columns, respectively.
- `x::Vector{Float64}`: The Raman shift values in cm⁻¹.
- `y::Vector{Float64}`: The intensity values corresponding to `x`.
- `multiple_spectra::Vector{<:Matrix{Float64}}`: A collection of spectra where each spectrum is a matrix with two columns (x and y values).
- `temperature_C::Float64`: Temperature in degrees Celsius.
- `laser_wavelength::Float64`: Wavelength of the excitation laser line in nanometers (nm).

# Options
- `correction::String="long"`: The equation used for the correction. Choose from:
  - `"long"`: Corrects using the Galeener equation with an additional ( nu_0^3 ) scaling factor for adimensionality (default).
  - `"galeener"`: Uses the original Galeener equation without the ( nu_0^3 ) factor.
  - `"hehlen"`: Applies the Hehlen equation, which includes density corrections to preserve low-frequency signals.
- `normalisation::String="area"`: Specifies whether to normalize the corrected signal. Options are:
  - `"area"`: Normalizes by integrating over the area under the curve.
  - `"intensity"`: Normalizes by peak intensity.
  - `"minmax"`: Scales between minimum and maximum values.
  - `"no"`: No normalization (default is `"area"`).
- `density::Float64=2210.0`: Density of the studied material in kg/m³, used only with the "hehlen" equation. Default is the density of silica glass.

# Returns
For single spectrum input (`x` and `y` or `spectrum`):
- `x_out::Vector{Float64}`: The Raman shift values (same as input x).
- `y_corr::Vector{Float64}`: The corrected intensity values.
- `ese_corr::Vector{Float64}`: The propagated errors on the corrected intensities.

For multiple spectra input (`multiple_spectra`):
- A vector of tuples where each tuple contains (`x_out`, `y_corr`, `ese_corr`) for a single spectrum.

# Notes
- The old API is not indicated but still available for backward compatibility: you can call tlcorrection(spectrum, ...) with spectrum = [x y]
- This correction uses the formula reported in Galeener and Sen (1978), Mysen et al. (1982), Brooker et al. (1988) and Hehlen et al. (2010).
- The "galeener" equation is the exact one reported in Galeener and Sen (1978), which is a modification from Shuker and Gammon (1970) for accounting of (vo - v)^4 dependence of the Raman intensity. See also Brooker et al. (1988) for further discussion.
- The "long" equation is that of Galeener and Sen (1978) corrected by a vo^3 coefficient for removing the cubic meter dimension of the equation of "galeener". This equation has been used in Mysen et al. (1982) and Le Losq et al. (2012).
- The "hehlen" equation is that reported in Hehlen et al. (2010). It actually originates before this publication (see Brooker et al. 1988). It uses a different correction that avoid crushing the signal below 500 cm-1. THerefore, it has the advantage of keeping intact the Boson peak signal in glasses.
# Notes
1. **Equations**:
   - The "long" equation is a modified version of Galeener's equation that includes a nu_0^3 scaling factor to remove cubic meter dimensions. This version has been widely used in studies such as Mysen et al. (1982) and Le Losq et al. (2012).
   - The "galeener" equation is the original form reported by Galeener and Sen (1978), which modifies Shuker and Gammon's (1970) approach to account for (vo - v)^4  dependence.
   - The "hehlen" equation, introduced by Hehlen et al. (2010), avoids signal suppression below 500 cm⁻¹, preserving features like the Boson peak in glasses.
2. **Error Propagation**:
   - Errors are calculated as sqrt{y} for raw data and propagated through the correction process.
3. **Normalization**:
   - If normalization is enabled, it is applied after the correction step

# Examples

## Example 1: Correcting a single spectrum, provided as x-y
```julia
using Spectra, Random
x = collect(100.0:1.0:1000)
y = rand(length(x)) # Example spectrum
temperature_C = 25.0
laser_wavelength = 532.0 # nm
x_out, y_corr, ese_corr = tlcorrection(x, y, temperature_C, laser_wavelength)
```
## Example 2: Correcting multiple spectra
```julia
spectrum1 = hcat(collect(100.0:10:1000), rand(91))
spectrum2 = hcat(collect(100.0:10:1000), rand(91))
multiple_spectra = [spectrum1, spectrum2]
corrected_spectra = tlcorrection(multiple_spectra, temperature_C, laser_wavelength)
```
# Errors
- Throws an error if an unsupported correction method is specified.
- Throws an error if normalization is not one of "area", "intensity", "minmax", or "no".
- Throws an error if input spectra have fewer than two columns or if no spectra are provided for multiple-spectra input.

# References

- Brooker et al. (1988) Journal of Raman Spectroscopy 19(2), 71-78.
- Galeener and Sen (1978) Physical Review B 17 (4): 1928–33.
- Hehlen (2010) Journal of Physics: Condensed Matter 22 (2): 025401.
- Le Losq et al. (2012) American Mineralogist, 97, 779–790.
- Mysen et al. (1982) American Mineralogist 67: 686–95.
- Shuker and Gammon (1970) Physical Review Letters 25 (4): 222–25.

"""
function tlcorrection(x::Vector{Float64}, y::Vector{Float64}, temperature_C::Float64, laser_wavelength::Float64;
	correction::String="long", normalisation::String="area", density::Float64=2210.0)

    h::Float64 = 6.626070040e-34   # J S    Plank constant from NIST
	hb::Float64 = 1.054571800e-34 # J S    Reduced Plank constant from NIST
    k::Float64 = 1.38064852e-23      # J K-1    Boltzman constant from NIST
    c::Float64 = 299792458.0         # M S-1    Speed of light from NIST
    nu0::Float64 = 1.0 ./laser_wavelength .*1.0e9     # nu0 laser is in M-1 (laser_wavelength is in nm)
    T::Float64 = temperature_C + 273.15    # K temperature
	# density is in KG M-3

    # Calculate the error on spectrum as sqrt(y). If y <= 0, then error = abs(y).
    ese = sqrt.(abs.(y))./abs.(y) # relative errors

	# get the Raman shift in m-1
	nu = 100.0.*x # cm-1 -> m-1 Raman shift

    # then we proceed to the correction
	if correction == "long"
		# Formula used in Mysen et al. (1982), Neuville and Mysen (1996) and Le Losq et al. (2012) 
		# (corrected for using the Planck constant in the last reference)
		# It is that reported in Brooker et al. (1988) with the addition of a scaling nu0^3 coefficient for adimentionality
    	frequency = nu0.^3.0.*nu./((nu0.-nu).^4) # frequency correction; dimensionless
    	boltzman = 1.0 .- exp.(-h.*c.*nu./(k.*T)) # temperature correction with Boltzman distribution; dimensionless
    	ycorr = y.*frequency.*boltzman; # correction

	elseif correction == "galeener"
		# This uses the formula reported in Galeener and Sen (1978) and Brooker et al. (1988); 
		# it uses the Bose-Einstein / Boltzman distribution
		# without the scaling nu0^3 coefficient
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

	if normalisation == "area" || normalisation == "intensity" || normalisation == "minmax"
		out = normalise(ycorr, x=x, method=normalisation)
	elseif normalisation == "no"
		out = copy(ycorr)
	else
		error("Set the optional normalisation parameter to area, intensity, minmax or no.")
	end

	# final error calculation
    ese_out = ese.*out 

    return x, out, ese_out
end
function tlcorrection(spectrum::Matrix{Float64}, temperature_C::Float64, laser_wavelength::Float64; correction::String="long", normalisation::String="area", density::Float64=2210.0)
	return tlcorrection(spectrum[:,1], spectrum[:,2], temperature_C, laser_wavelength, correction=correction, normalisation=normalisation, density=density)
end
function tlcorrection(multiple_spectra::Vector{<:Matrix{Float64}}, temperature_C::Float64, laser_wavelength::Float64; correction::String="long", normalisation::String="area", density::Float64=2210.0)

	# check if the input is a vector of spectra
	if length(multiple_spectra) == 0
		error("No spectra provided.")
	end

	# apply tlcorrection to each spectrum in the vector
	corrected_spectra = [tlcorrection(spectrum, temperature_C, laser_wavelength, correction=correction, normalisation=normalisation, density=density) for spectrum in multiple_spectra]

	return corrected_spectra
end