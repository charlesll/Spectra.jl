"""
    fit_peaks(x, y, peaks_info; sigma=nothing, backend=:qGN, relax=4, maxiter=100)

Fit a sum of peaks (Gaussian, Lorentzian, or Pseudo-Voigt) to observed data `y` over the domain `x`. 

This function allows for constrained optimization within specific parameter bounds and incorporates prior uncertainties into the fitting process. The misfit function ``S`` used is (Tarantola, 2005, eq. 3.46):

`` S = \\sum(\\frac{(y - y_{calc})^2}{\\sigma_y^2}) + \\frac{(m_{current} - m_{prior})^2}{\\sigma_m_{prior}^2}) ``

with ``y_{calc}`` the calculated signal, ``y`` the observed signal, ``m_{current}`` the current model parameters, and ``m_{prior}`` the prior model parameters. The first term is the data misfit, and the second term is the model prior misfit.

Two optimization backends are available:

1. **`:qGN`: Quasi Gauss-Newton backend (recommended)**:
   - Implements the quasi Gauss-Newton algorithm described in Tarantola (2005, eq. 3.51).
   - Includes a relaxation parameter `relax` for controlling step sizes.
   - Fast.

2. **`:Optim`: Optim.jl backend**:
   - BFGS algorithm with automatic differentiation.
   - Less fast...
   - ...but uses FMinxBox so you can further constrain the fit within specific boundaries.

# Arguments

- `x::Vector{Float64}`: The independent variable (domain) over which the peaks are defined.
- `y::Vector{Float64}`: The observed data to be fitted.
- `peaks_info::Vector{Tuple{Symbol, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}}}`:
  A list of tuples specifying the peak types and their parameters:
  - `Symbol`: Type of peak (`:gaussian`, `:lorentzian`, `:pseudovoigt`, `:pearson7`).
  - `Vector{Float64}`: Initial parameter values:
    - for `:gaussian` or `lorentzian`: [amplitude, center, width]. 
    - for `:pseudovoigt`: [amplitude, center, width, fraction].
    - for `:pearson7`: [amplitude, center, width, exponent].
  - `Vector{Float64}`: Prior uncertainties for each parameter.
  - `Vector{Float64}`: Lower bounds for each parameter (warning: only used when `backend=:Optim`).
  - `Vector{Float64}`: Upper bounds for each parameter (warning: only used when `backend=:Optim`).

# Keyword Arguments

- `sigma::Union{Nothing, Vector{Float64}}`: Optional uncertainties on the observed data (`y`). If provided, these are used to weight the misfit function.
- `backend::Symbol`: Specifies the optimization backend (`:Optim` or `:qGN`).
- `relax::Int`: Relaxation parameter for the quasi Gauss-Newton algorithm (default is 4). The larger it is, the smaller the steps, the better the numerical stability, but the longer the fit. If the fit fails, increase this value.
- `maxiter::Int`: Maximum number of iterations for the quasi Gauss-Newton algorithm (default is 100). If you increase `relax`, consider also increasing `maxiter`.

# Returns

A named tuple containing:
- `fit`: The optimization result object from Optim.jl or the final fitted parameters from the quasi Gauss-Newton algorithm.
- `peak_results`: A list of fitted peak parameters, their uncertainties and their areas.
- `peak_areas`: A vector containing the peak area fractions.
- `residuals`: The difference between the observed data (`y`) and the fitted model predictions.
- `hessian`: The Hessian matrix at the optimum point.
- `CMpost`: The model parameters covariance matrix at the optimum point.
- `plot_fit`: A function to plot the observed data and fitted model.
- `print_params`: A function to print fitted parameters and uncertainties.

# Note
- There is a fine balance between `relax` and `maxiter` for the quasi-Newton method. If you want speed, try setting `relax` to the smallest possible value, like 1-2 and `maxiter` to 20. Check the convergence by comparing results with runs with higher `relax` and `maxiter` (such as 100 and 1000, respectively).
- With the quasi-Newton method, if you run into numerical errors, like complex numbers and so on, increase both `relax` and `maxiter`. Also check your errors on your signal. They should be realistic.

# Examples

We fit a signal composed of three Gaussian, Lorentzian and Pseudovoigt peaks.
The params_uncertainties vector constrain the 
fit around the provided prior values. You can set broad uncertainties if you don't want to constrain the fit.

## Example 1: quasi-Newton method

```@example
x_fit = collect(0:0.2:100)
y_fit_perfect = (gaussian(x_fit, [10.0, 35.0, 10.0]) 
+ lorentzian(x_fit, [15.0, 55.0, 3.0]) 
+ pseudovoigt(x_fit, [20.0, 45.0, 2.0, 0.4]))
noise = randn(length(x_fit))*0.3
y_fit = y_fit_perfect + noise

peaks_info = [
    # (peak_type, prior_params, prior_uncertainties, lower_bounds, upper_bounds)
    (:gaussian, [10.5, 33.0, 10.0], [5.0, 10.0, 5.0], [0., 0., 0.], [Inf, Inf, 50.]),
    (:lorentzian, [15.5, 57.0, 3.1], [5.0, 10., 3.0], [0., 0., 0.], [100., 100., 50.]),
    (:pseudovoigt, [16.5, 44.0, 3.0, 0.4], [10.0, 10.0, 3.0, 0.3], [0., 0., 0., 0.], [100., 100., 50., 1.0])
]

result = fit_peaks(x_fit, y_fit, peaks_info, sigma=noise, backend=:qGN, relax=3, maxiter=100)
result.print_params()
```

If needed, you can constrain further the search by decreasing the uncertainties on
prior values for parameters that you think you know well. 
Let's do that for the pseudovoigt `lorentzian_fraction` for the sack of example:
```@example
peaks_info = [
# (peak_type, initial_params, uncertainties, lower_bounds, upper_bounds)
(:gaussian, [10.5, 30.0, 3.0], [5.0, 10.0, 5.0], [0., 25., 0.], [Inf, Inf, 50.]),
(:lorentzian, [20.5, 55.0, 3.0], [10.0, 10., 5.0], [0., 50., 0.], [100., 100., 50.]),
(:pseudovoigt, [20.5, 44.0, 3.0, 0.4], [5.0, 5., 5.0, 0.01], [0., 0., 0., 0.], [100., 100., 50., 1.])
]

result = fit_peaks(x_fit, y_fit, peaks_info, sigma=noise, backend=:qGN, relax=3, maxiter=100)
result.print_params()
```

## Example 2: using Optim.jl backend

The `lower_bounds` and `upper_bounds` are only used by Optim.jl to constrain the fit within boundaries. Here for instance the `lorentzian_fraction` is comprised between 0 and 1, or the 
peak hwhms are constraind to be positive. To use them, use the `:Optim` backend:
```@example
result = fit_peaks(x_fit, y_fit, peaks_info, sigma=noise, backend=:Optim)
result.print_params(digits=3)
````
## Example 3: Plot the fit

After a fit, you can plot the signal, fit and components using:
```@example
result.plot_fit(components=true)
```

# References

Tarantola A., *Inverse Problem Theory and Methods for Model Parameter Estimation*, SIAM (2005), Chapter 3.

"""
function fit_peaks(x, y, peaks_info; 
    sigma=nothing, 
    backend=:Optim, 
    relax=4, 
    maxiter=100)

    param_indices = Dict()
    all_params_prior = Float64[]
    all_prior_uncertainties = Float64[]
    all_lower_bounds = Float64[]
    all_upper_bounds = Float64[]
    current_idx = 1

    for (i, (peak_type, params, uncertainties, lower_bounds, upper_bounds)) in
        enumerate(peaks_info)
        n_params = peak_type == :pseudovoigt || peak_type == :pearson7 ? 4 : 3
        param_indices[i] = current_idx:(current_idx + n_params - 1)
        append!(all_params_prior, params)
        append!(all_prior_uncertainties, uncertainties)
        append!(all_lower_bounds, lower_bounds)
        append!(all_upper_bounds, upper_bounds)
        current_idx += n_params
    end

    # Define peak function dispatch dictionary
    peak_functions = Dict(
        :gaussian => gaussian,
        :lorentzian => lorentzian,
        :pseudovoigt => pseudovoigt,
        :pearson7 => pearson7,
    )

    # need for speed : if sigma is not provided, we set it to an array of 1 to avoid
    # else-if terms in the functions that participate in the optimization
    if sigma == nothing
        sigma = ones(size(y,1))
    end

    # define the covariance matrices in advance
    CM = diagm(all_prior_uncertainties .^ 2)
    CD = diagm(sigma .^ 2)
    ICM = inv(CM) # and its inverse
    ICD = inv(CD) # and its inverse

    # define in advance the loss distribution used in :Optim
    data_loss_dist = MvNormal(y, CD)
    mprior_loss_dist = MvNormal(all_params_prior, CM)

    # Direct model
    function model(params)
        result = zero(x) * zero(params[1])  # creates a vector of zeros with the right type
        for (i, (peak_type, _, _, _, _)) in enumerate(peaks_info)
            p = params[param_indices[i]]
            result .+= peak_functions[peak_type](x, p)
        end
        return result
    end
    function model_lsqfit(x, params)
        return model(params)
    end

    """
        optimize_GN(params; maxiter=maxiter, relax=relax)
    
    quasi-Newton optimizer, equation 3.51, page 69, Tarantola 2005
    """
    function optimize_GN(params; maxiter=maxiter, relax=relax)
        
        #the model à priori, as a vector
        mprior = deepcopy(params)
        mcurrent = deepcopy(params) # initialize the algorithm at the à priori model values

        for i = 1:maxiter
            # data and model residuals
            dataresiduals = model(mcurrent) - y
            modelresiduals = mcurrent-mprior

            # derivatives
            G = ForwardDiff.jacobian(model, mcurrent) 
            GT = transpose(G)

            # equation 3.51, page 69, Tarantola 2005
            gradient = GT*ICD*dataresiduals + CM\modelresiduals # The hessian gradient
            direction = (GT*ICD*G+ICM)\gradient # The direction
            mcurrent = mcurrent - direction ./ relax # update mcurrent with a relaxation parameter

        end

        # derivatives at the optimum
        G = ForwardDiff.jacobian(model, mcurrent)
        GT = transpose(G)

        # formula 3.53, Tarantola 2005
        CMpost = CM - CM*GT*inv(G*CM*GT + CD)*G*CM
        errors = sqrt.(diag(CMpost))

        return mcurrent, errors, CMpost
    end

    # loss function with priors:
    # log likelihood + log prior
    function loss(params)
        data_loss = -logpdf(data_loss_dist,  model(params))
        mprior_loss = -logpdf(mprior_loss_dist, params)
        return data_loss + mprior_loss
    end
    function loss_details(params) # a dev version to get access to the data/prior losses
        data_loss = -logpdf(data_loss_dist,  model(params))
        mprior_loss = -logpdf(mprior_loss_dist, params)
        return data_loss, mprior_loss
    end
    # loss by direct least-squares calculation
    # does not affect convergence, but do not use for Hessian calculation
    function ls_loss(params)
        data_loss = sum((model(params).-y).^2 ./ sigma.^2)
        mprior_loss = sum((all_params_prior.-params).^2 ./ all_prior_uncertainties.^2)
        out = data_loss + mprior_loss
        return out
    end

    # Perform the optimization
    if backend == :Optim
        result = optimize(
            ls_loss, # we use the ls_loss => quicker, does not affect convergence
            all_lower_bounds,
            all_upper_bounds,
            all_params_prior,
            #LBFGS();
            Fminbox(LBFGS());
            autodiff = :forward
        )
        fitted_params = Optim.minimizer(result)
        # Calculate errors (approximate)
        # We use now the loss() function that result in mathematically-sound Hessian calculation
        # because of the scaling involved in the log likelihood and prior
        hessian_matrix = ForwardDiff.hessian(loss, fitted_params)
        CMpost = inv(hessian_matrix)
        errors = sqrt.(diag(CMpost))
    elseif backend == :qGN
        fitted_params, errors, CMpost = optimize_GN(
            all_params_prior; maxiter=maxiter, relax=relax
        )
        hessian_matrix = ForwardDiff.hessian(loss, fitted_params)
    elseif backend == :LsqFit
        fit = curve_fit(model_lsqfit, x, y, all_params_prior)
        hessian_matrix = nothing
        CMpost = nothing
        fitted_params = fit.param
        errors = stderror(fit)
    else
        error("Backend not implemented")
    end

    residuals = y - model(fitted_params)

    # Organize results by peak
    peak_results = []
    normalised_peak_areas = Float64[]
    for (i, (peak_type, _, _)) in enumerate(peaks_info)
        # get the parameters for the peak we are interested in
        params_ = fitted_params[param_indices[i]]
        amplitude_ = params_[1]
        hwhm_ = params_[3]
        lorentzian_fraction_ = peak_type == :pseudovoigt ? params_[4] : nothing
        exponent_ = peak_type == :pearson7 ? params_[4] : nothing

        # calculate the area and store it in a vector of the normalised peak areas
        area_peak_out = area_peaks(
            peak_type;
            amplitude=amplitude_,
            hwhm=hwhm_,
            lorentzian_fraction=lorentzian_fraction_,
            exponent=exponent_,
        )
        push!(normalised_peak_areas, area_peak_out)

        # store the peak results also in peak_results Named Tuple
        push!(
            peak_results,
            (
                peak_type=peak_type,
                params=params_,
                errors=errors[param_indices[i]],
                area=area_peak_out,
            ),
        )
    end

    """
        print_params(digits=4)
                
    print the parameters after the fit

    """
    function print_params(;digits=4)
        for (i, peak) in enumerate(peak_results)
            println("Peak $i ($(peak.peak_type)):")
            if peak.peak_type == :pseudovoigt
                param_names = ["amplitude", "center", "width", "fraction"]
            elseif peak.peak_type == :pearson7
                param_names = ["amplitude", "center", "width", "exponent"]
            else
                param_names = ["amplitude", "center", "width"]
            end

            for j in 1:length(peak.params)
                println(
                    "  $(param_names[j]): $(round(peak.params[j], digits=digits)) ± $(round(peak.errors[j], digits=digits))",
                )
            end
        end
    end
    """
        plot_fit(; components=false)

    plot the data and the fit

    """
    function plot_fit(; components=false)
        p = scatter(x, y; label="Data", legend=:topleft)
        plot!(p, x, model(fitted_params); lw=2, label="Fit")

        if components
            for (i, peak) in enumerate(peak_results)
                component_params = fitted_params[param_indices[i]]
                if peak.peak_type == :gaussian
                    plot!(p, x, gaussian(x, component_params); label="Gaussian $i")
                elseif peak.peak_type == :lorentzian
                    plot!(p, x, lorentzian(x, component_params); label="Lorentzian $i")
                elseif peak.peak_type == :pseudovoigt
                    plot!(p, x, pseudovoigt(x, component_params); label="PseudoVoigt $i")
                elseif peak.peak_type == :pearson7
                    plot!(p, x, pearson7(x, component_params); label="Pearson7 $i")
                end
            end
        end

        return p
    end

    # Return fit results and plotting function
    return (
        fit=result,
        peak_results=peak_results,
        peak_areas=normalised_peak_areas/sum(normalised_peak_areas),
        residuals=residuals,
        hessian=hessian_matrix,
        CMpost=CMpost,
        plot_fit=plot_fit,
        print_params=print_params,
        loss_details=loss_details,
        fitted_params = fitted_params,
        param_indices = param_indices
    )
end

"""
    print_params(peak_results; digits=4)
            
    print the parameters after the fit
"""
function print_params(peak_results; digits=4)
    for (i, peak) in enumerate(peak_results)
        println("Peak $i ($(peak.peak_type)):")
        if peak.peak_type == :pseudovoigt
            param_names = ["amplitude", "center", "width", "fraction"]
        elseif peak.peak_type == :pearson7
            param_names = ["amplitude", "center", "width", "exponent"]
        else
            param_names = ["amplitude", "center", "width"]
        end

        for j in 1:length(peak.params)
            println(
                "  $(param_names[j]): $(round(peak.params[j], digits=digits)) ± $(round(peak.errors[j], digits=digits))",
            )
        end
    end
end