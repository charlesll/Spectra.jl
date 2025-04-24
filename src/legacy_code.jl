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
function xshift_correction(
    full_x::Array{Float64},
    full_shifted_y::Array{Float64},
    ref_x::Array{Float64},
    ref_y::Array{Float64},
    shifted_y::Array{Float64},
)
    fit = curve_fit(xshift_inversion, [ref_x shifted_y], ref_y, [1.0, 1.0, 1.0])
    parameters = fit.param
    return xshift_direct(full_x, full_shifted_y, parameters[1])
end

"""
	xshift_inversion(data::Array{Float64},p::Array{Float64})

for xshift_correction to correct for shifts in X between two spectra

"""
function xshift_inversion(data::Array{Float64}, p::Array{Float64})
    xaxis = data[:, 1]
    shifted1 = data[:, 2]
    spl = Spline1D(xaxis-p[1], shifted1 .* p[2] + shifted1 .^ 2.0 .* p[3])
    y = spl(xaxis)
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
function pearson7(
    a1::Array{Float64},
    a2::Array{Float64},
    a3::Array{Float64},
    a4::Array{Float64},
    x::Array{Float64};
    style::String="None",
)
    segments = zeros(size(x)[1], size(a1)[1])
    if style == "None"
        for i in 1:size(a1)[1]
            segments[:, i] =
                a1[i] ./ (
                    1.0 .+
                    ((x[:, 1] .- a2[i]) ./ a3[i]) .^ 2 .* (2.0 .^ (1.0 ./ a4[i]) .- 1.0)
                )
        end
    elseif style == "poly"
        for i in 1:size(a1)[1]
            segments[:, i] =
                poly(vec(a1[i, :]), x[:, 2]) ./ (
                    1.0 .+
                    (
                        (x[:, 1] .- (poly(vec(a2[i, :]), x[:, 2]))) ./
                        poly(vec(a3[i, :]), x[:, 2])
                    ) .^ 2 .* (2.0 .^ (1, 0 ./ vec(a4[i, :])) .- 1.0)
                )
        end
    else
        error("Not implemented, see documentation")
    end
    return sum(segments; dims=2), segments
end

"""

	bandarea(Amplitude::Array{Float64},HWHM::Array{Float64}; peak_shape = "Gaussian", error_switch = "no", eseAmplitude::Array{Float64} = [0.0], eseHWHM::Array{Float64} = [0.0])

This function allows calculating the area under a specific band, with different shapes. For now, only Gaussian bands are supported, but other band shapes will be added soon.

Inputs
------

	Amplitude: Array{Float64}
		peak amplitude
	HWHM: Array{Float64}
		peak half width at half maximum

Options
-------

	peak_shape: String
		shape of the peak. Only "Gaussian" is supported for now
	error_switch: String
		should be "yes" or "no". If "yes", the arrays containing the errors affecting the band amplitude and widhts should be provided in eseAmplitude and eseHWHM (see below);
	eseAmplitude: Array{Float64}
		array containing the errors affecting Amplitude
	eseHWHM: Array{Float64}
		array containing the errors affecting HWHM;

Outputs
-------
	area: Array{Float64}
		array containing peak areas

	if error_switch is set to "yes", then a second output is provided:

	esearea: Array{Float64}
		array that contains the propagated errors affecting the areas calculations.
"""
function bandarea(
    Amplitude::Array{Float64},
    HWHM::Array{Float64};
    peak_shape="Gaussian",
    error_switch="no",
    eseAmplitude::Array{Float64}=[0.0],
    eseHWHM::Array{Float64}=[0.0],
)

    # first we check the desired peak shape, and apply the relevant calculation
    if peak_shape == "Gaussian"
        area::Array{Float64} = sqrt(pi ./ log(2)) .* Amplitude .* HWHM # Gaussian area, HWHM is the half-width at half-maximum
        if error_switch == "yes" # if errors are desired, perform the error calculation
            if size(eseAmplitude) != size(Amplitude) || size(eseHWHM) != size(HWHM) # error check
                error(
                    "Please check that you entered good arrays for the errors on the amplitude and widths of the bands",
                )
            end
            esearea::Array{Float64} = sqrt(
                (pi ./ log(2) .* HWHM) .^ 2 .* eseAmplitude .^ 2 +
                (pi ./ log(2) .* Amplitude) .^ 2 .* eseHWHM .^ 2,
            )
        end
    else
        error("Not yet implemented.")
    end

    # Depending on the error switch, we output only the areas or also the errors
    if error_switch == "no"
        # a quick test to see if the output is an array or should be a Float64 number
        if size(area) == (1,)
            return area[1] # we return a Float64 number
        else
            return area # we return an array
        end
    elseif error_switch == "yes"
        # a quick test to see if the output is an array or should be 2 Float64 numbers
        if size(area) == (1,)
            return area[1], esearea[1] # we return two Float64 numbers
        else
            return area, esearea # we return an array
        end
    else
        error(
            "error_switch should be equal to yes or no. Please select the appropriate value.",
        )
    end
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
function gaussiennes(
    amplitude::Array{Float64},
    centre::Array{Float64},
    hwhm::Array{Float64},
    x::Array{Float64};
    style::String="None",
)
    segments = zeros(size(x)[1], size(amplitude)[1])
    if style == "None"
        for i in 1:size(amplitude)[1]
            segments[:, i] =
                amplitude[i] .* exp.(-log.(2.0) .* ((x[:, 1] .- centre[i]) ./ hwhm[i]) .^ 2)
        end
    elseif style == "poly"
        for i in 1:size(amplitude)[1]
            segments[:, i] =
                poly(vec(amplitude[i, :]), x[:, 2]) .* exp.(
                    -log.(2) .*
                    (
                        (x[:, 1] .- (poly(vec(centre[i, :]), x[:, 2]))) ./
                        poly(vec(hwhm[i, :]), x[:, 2])
                    ) .^ 2,
                )
        end
    else
        error("Not implemented, see documentation")
    end
    return sum(segments; dims=2), segments
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
function lorentziennes(
    amplitude::Array{Float64},
    centre::Array{Float64},
    hwhm::Array{Float64},
    x::Array{Float64};
    style::String="None",
)
    segments = zeros(size(x)[1], size(amplitude)[1])
    if style == "None"
        for i in 1:size(amplitude)[1]
            segments[:, i] =
                amplitude[i] ./ (1.0 .+ ((x[:, 1] .- centre[i]) ./ hwhm[i]) .^ 2)
        end
    elseif style == "poly"
        for i in 1:size(amplitude)[1]
            segments[:, i] =
                poly(vec(amplitude[i, :]), x[:, 2]) ./ (
                    1.0 .+
                    (
                        (x[:, 1] .- (poly(vec(centre[i, :]), x[:, 2]))) ./
                        poly(vec(hwhm[i, :]), x[:, 2])
                    ) .^ 2
                )
        end
    else
        error("Not implemented, see documentation")
    end
    return sum(segments; dims=2), segments
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
function pseudovoigts(
    amplitude::Array{Float64},
    centre::Array{Float64},
    hwhm::Array{Float64},
    lorentzian_fraction::Array{Float64},
    x::Array{Float64};
    style::String="None",
)
    segments = zeros(size(x)[1], size(amplitude)[1])
    test1 = findall(lorentzian_fraction .< 0.0)
    test2 = findall(lorentzian_fraction .> 1.0)
    if size(test1)[1] != 0 || size(test2)[1] != 0
        error("lorentzian_fraction should be comprised between 0 and 1")
    end
    if style == "None"
        for i in 1:size(amplitude)[1]
            segments[:, i] =
                amplitude[i] .* (
                    (1.0 .- lorentzian_fraction[i]) .*
                    exp.(-log(2) .* ((x[:, 1] .- centre[i]) ./ hwhm[i]) .^ 2) .+
                    lorentzian_fraction[i] ./
                    (1.0 .+ ((x[:, 1] .- centre[i]) ./ hwhm[i]) .^ 2)
                )
        end
    elseif style == "poly"
        for i in 1:size(amplitude)[1]
            segments[:, i] =
                vec(lorentzian_fraction[i, :]) .* (
                    poly(vec(amplitude[i, :]), x[:, 2]) ./ (
                        1 .+
                        (
                            (x[:, 1] .- (poly(vec(centre[i, :]), x[:, 2]))) ./
                            poly(vec(hwhm[i, :]), x[:, 2])
                        ) .^ 2
                    )
                ) .+
                (1.0 .- vec(lorentzian_fraction[i, :])) .*
                poly(vec(amplitude[i, :]), x[:, 2]) .* exp.(
                    -log(2) .*
                    (
                        (x[:, 1] .- (poly(vec(centre[i, :]), x[:, 2]))) ./
                        poly(vec(hwhm[i, :]), x[:, 2])
                    ) .^ 2,
                )
        end
    else
        error("Not implemented, see documentation")
    end
    return sum(segments; dims=2), segments
end

"""
	bootperf(params_boot::Array{Float64}; plotting::String = "True", parameter::Int64 = 1, feature::Int64 = 1, histogram_step::Int64 = 100, savefigures::String = "False", save_bootrecord::String = "Boot_record.pdf", save_histogram::String = "Boot_histogram.pdf")

Inputs
------

	params_boot[i,j] or [i,j,k]
		array with i the bootstrrap experiment, j the parameter and k the feature being calculated (e.g., a peak during peak fitting)
	plotting: String
		switch to plotting mode ("True") or not ("False"). If true, parameter and feature must be provided, otherwise an error message is returned.

	histogram_step
		integer value to control the histogram X axis division.
	savefigures: "True" or "False"
		explicit, save in the current working directory.
	save_bootrecord:
		Name for the graphic showing the bootstrap mean and std evolutions, with extension.
	save_histogram:
		Name for the graphic showing the histogram for the parameter and feature of interest, with extension.

Outputs
-------

	std_record, mean_record, the arrays recording how the standard deviation and mean of the parameters as a function of the bootstrap advance.

"""
function bootperf(
    params_boot::Array{Float64};
    plotting::String="True",
    parameter::Int64=1,
    feature::Int64=1,
    histogram_step::Int64=100,
    savefigures::String="False",
    save_bootrecord::String="Boot_record.pdf",
    save_histogram::String="Boot_histogram.pdf",
)

    # test
    if plotting != "True" && plotting != "False"
        error("Error: plotting should be set to True or False")
    end

    nb_boot::Int64 = size(params_boot)[1]
    nb_param::Int64 = size(params_boot)[2]

    if ndims(params_boot) == 2 #### TWO DIMENSIONS PARAMS_BOOT
        mean_record::Array{Float64} = zeros(size(params_boot))
        std_record::Array{Float64} = zeros(size(params_boot))
        for k in 1:nb_boot
            for j in 1:size(params_boot)[2]
                if k == 1
                    mean_record[k, j] = params_boot[k, j]
                    std_record[k, j] = 0
                else
                    mean_record[k, j] = sum(mean(params_boot[1:k, j], 1))
                    std_record[k, j] = sum(std(params_boot[1:k, j], 1))
                end
            end
        end

        if plotting == "True" # Graphics
            # tests:
            if savefigures != "True" && savefigures != "False"
                error("savefigure should be set to True or False.")
            end
            if isa(save_bootrecord, String) != true
                error(
                    "Error: save_bootrecord is not a valid savename. Check it is a string of characters.",
                )
            end
            if isa(save_histogram, String) != true
                error(
                    "Error: save_histogram is not a valid savename. Check it is a string of characters.",
                )
            end

            if parameter == 0 || histogram_step == 0
                error("parameter or histogram_step appear to have been set to 0...")
            elseif parameter > size(params_boot)[2]
                error("Value entered for parameter is too high.")
            end

            fig = figure()
            #title("Mean and Std of parameter $(parameter), feature $(feature)",fontsize = 10, fontweight = "bold", fontname="Arial")
            subplot(121)
            plot(1:nb_boot, mean_record[:, parameter]; label="Mean")
            xlabel(
                "Number of iterations during bootstrap";
                fontsize=10,
                fontweight="bold",
                fontname="Arial",
            )
            ylabel("Mean value"; fontsize=10, fontweight="bold", fontname="Arial")

            subplot(122)
            plot(1:nb_boot, std_record[:, parameter]; label="Std")
            xlabel(
                "Number of iterations during bootstrap";
                fontsize=10,
                fontweight="bold",
                fontname="Arial",
            )
            ylabel("Standard deviation"; fontsize=10, fontweight="bold", fontname="Arial")

            if savefigures == "True"
                savefig(save_bootrecord)
            end

            fig2 = figure()
            h = plt[:hist](params_boot[:, parameter], 100)
            #title("Histogram of parameter $(parameter), feature $(feature)",fontsize = 10, fontweight = "bold", fontname="Arial")
            xlabel("Parameter value"; fontsize=10, fontweight="bold", fontname="Arial")
            ylabel(
                "Number of bootstrap values";
                fontsize=10,
                fontweight="bold",
                fontname="Arial",
            )

            if savefigures == "True"
                savefig(save_histogram)
            end
        end # end if plotting

        return mean_record, std_record

    elseif ndims(params_boot) == 3 #### THREE DIMENSIONS PARAMS_BOOT
        mean_record = zeros(size(params_boot))
        std_record = zeros(size(params_boot))
        for k in 1:nb_boot
            for j in 1:size(params_boot)[2]
                for i in 1:size(params_boot)[3]
                    if k == 1
                        mean_record[k, j, i] = params_boot[k, j, i]
                        std_record[k, j, i] = 0
                    else
                        mean_record[k, j, i] = sum(mean(params_boot[1:k, j, i], 1))
                        std_record[k, j, i] = sum(std(params_boot[1:k, j, i], 1))
                    end
                end
            end
        end

        if plotting == "True" # Graphics
            # tests:
            if savefigures != "True" && savefigures != "False"
                error("savefigure should be set to True or False.")
            end
            if isa(save_bootrecord, String) != true
                error(
                    "Error: save_bootrecord is not a valid savename. Check it is a string of characters.",
                )
            end
            if isa(save_bootrecord, String) != true
                error(
                    "Error: save_histogram is not a valid savename. Check it is a string of characters.",
                )
            end

            if parameter == 0 || feature == 0 || histogram_step == 0
                error(
                    "parameter, feature or histogram_step appear to have been set to 0..."
                )
            elseif parameter > size(params_boot)[2]
                error("Value entered for parameter is too high.")
            elseif feature > size(params_boot)[3]
                error("Value entered for feature is too high.")
            end

            fig = figure()
            #title("Mean and Std of parameter $(parameter), feature $(feature)",fontsize = 10, fontweight = "bold", fontname="Arial")
            subplot(121)
            plot(1:nb_boot, mean_record[:, parameter, feature]; label="Mean")
            xlabel(
                "Number of iterations during bootstrap";
                fontsize=10,
                fontweight="bold",
                fontname="Arial",
            )
            ylabel("Mean value"; fontsize=10, fontweight="bold", fontname="Arial")

            subplot(122)
            plot(1:nb_boot, std_record[:, parameter, feature]; label="Std")
            xlabel(
                "Number of iterations during bootstrap";
                fontsize=10,
                fontweight="bold",
                fontname="Arial",
            )
            ylabel("Standard deviation"; fontsize=10, fontweight="bold", fontname="Arial")

            if savefigures == "True"
                savefig(save_bootrecord)
            end

            fig2 = figure()
            h = plt[:hist](params_boot[:, parameter, feature], 100)
            #title("Histogram of parameter $(parameter), feature $(feature)",fontsize = 10, fontweight = "bold", fontname="Arial")
            xlabel("Parameter value"; fontsize=10, fontweight="bold", fontname="Arial")
            ylabel(
                "Number of bootstrap values";
                fontsize=10,
                fontweight="bold",
                fontname="Arial",
            )

            if savefigures == "True"
                savefig(save_histogram)
            end
        end # end if plotting

        return mean_record, std_record

    elseif ndims(params_boot) > 3
        error("Not implemented, params_boot should have 3 dimensions, maximum.")
    end
end # end function
