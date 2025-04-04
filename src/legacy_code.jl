# Functions that were in Spectra but are not usefull or actually not finished... I remove them

"""
	xshift_correction(full_x::Array{Float64}, full_shifted_y::Array{Float64}, ref_x::Array{Float64}, ref_y::Array{Float64},shifted_y::Array{Float64})

To correct a shift between two spectra using a reference peak.

Inputs
------

	full_x: Array{Float64}
		x values that are not good
	full_shifted_y: Array{Float64}
		y values associated with full_x
	ref_x: Array{Float64}
		x values that are good
	ref_y: Array{Float64}
		y values associated with ref_x
	shifted_y: Array{Float64}
		y values associated with a selected range of full_x that corresponds to ref_x (for instance, a specific peak that you want to use to correct the shift).

Outputs
-------

	full_x: Array{Float64}
		same as input
	corrected_y: Array{Float64}
		the full_shifted_y values corrected from the shift
	p: Array{Float64}
		same as input.

ref_x is the common X axis of two particular ref_y and shifted_y signals, that should be for instance an intense and well defined peak in your spectra. If ref_y and shifted_y do not share the same X axis, you can use first the Dierckx spline to re-sample one of them and have both sharing a common X axis. See the examples for further details.
"""
function xshift_correction(full_x::Array{Float64}, full_shifted_y::Array{Float64}, ref_x::Array{Float64}, ref_y::Array{Float64},shifted_y::Array{Float64})
	fit = curve_fit(xshift_inversion, [ref_x shifted_y], ref_y, [1.0, 1.0,1.0])
	parameters = fit.param
	return xshift_direct(full_x, full_shifted_y, parameters[1])
end

"""
	xshift_inversion(data::Array{Float64},p::Array{Float64})

for xshift_correction to correct for shifts in X between two spectra

"""
function xshift_inversion(data::Array{Float64},p::Array{Float64})
    xaxis = data[:,1]
    shifted1 = data[:,2]
    spl = Spline1D(xaxis-p[1], shifted1.*p[2] + shifted1.^2.0.*p[3])
    y = spl(xaxis)
end