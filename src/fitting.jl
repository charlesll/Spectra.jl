"""
    fit_peaks(x, y, peaks_info; sigma=nothing, backend=:Optim, relax=100, maxiter=100)

Fit a sum of peaks (Gaussian, Lorentzian, or Pseudo-Voigt) to observed data `y` over the domain `x`. 

This function allows for constrained optimization within specific parameter bounds and incorporates prior uncertainties into the fitting process. The misfit function used is:

    misfit = sum(((predicted - y)^2) / (sigma^2)) +
             sum(((params - params_prior)^2) / (uncertainty_params_prior^2))

This corresponds to equation 3.46 in Tarantola (2005). Two optimization backends are available:

1. **Optim.jl backend**:
   - Uses the LBFGS algorithm with automatic differentiation.
   - Suitable for general-purpose optimization.

2. **Quasi Gauss-Newton backend**:
   - Implements the quasi Gauss-Newton algorithm described in Tarantola (2005, eq. 3.51).
   - Includes relaxation parameters for controlling step sizes.

# Arguments

- `x::Vector{Float64}`: The independent variable (domain) over which the peaks are defined.
- `y::Vector{Float64}`: The observed data to be fitted.
- `peaks_info::Vector{Tuple{Symbol, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}}}`:
  A list of tuples specifying the peak types and their parameters:
  - `Symbol`: Type of peak (`:gaussian`, `:lorentzian`, or `:pseudovoigt`).
  - `Vector{Float64}`: Initial parameter values ([amplitude, center, width] or [amplitude, center, width, fraction]).
  - `Vector{Float64}`: Prior uncertainties for each parameter.
  - `Vector{Float64}`: Lower bounds for each parameter.
  - `Vector{Float64}`: Upper bounds for each parameter.

# Keyword Arguments

- `sigma::Union{Nothing, Vector{Float64}}`: Optional uncertainties on the observed data (`y`). If provided, these are used to weight the misfit function.
- `backend::Symbol`: Specifies the optimization backend (`:Optim` or `:qGN`).
- `relax::Int`: Relaxation parameter for the quasi Gauss-Newton algorithm (default is 100).
- `maxiter::Int`: Maximum number of iterations for the quasi Gauss-Newton algorithm (default is 100).

# Returns

A named tuple containing:
- `fit`: The optimization result object from Optim.jl or the final fitted parameters from the quasi Gauss-Newton algorithm.
- `peak_results`: A list of fitted peak parameters and their uncertainties.
- `residuals`: The difference between the observed data (`y`) and the fitted model predictions.
- `hessian`: The Hessian matrix at the optimum point (if available).
- `plot_fit`: A function to plot the observed data and fitted model.
- `print_params`: A function to print fitted parameters and uncertainties.
- `CMpost`: Posterior covariance matrix for parameter uncertainties (only available with the quasi Gauss-Newton backend).

# Examples

```julia
x_fit = sort(rand(1000)*100)
y_fit_perfect = (gaussian(x_fit, [10.0, 35.0, 10.0]) 
+ lorentzian(x_fit, [15.0, 55.0, 3.0]) 
+ pseudovoigt(x_fit, [20.0, 45.0, 2.0, 0.4]))
noise = randn(length(x_fit))*0.3
y_fit = y_fit_perfect + noise

peaks_info = [
# (type, initial_params, uncertainties, lower_bounds, upper_bounds)
(:gaussian, [10.5, 30.0, 3.0], [5.0, 10.0, 5.0], [0., 0., 0.], [Inf, Inf, 50.]),
(:lorentzian, [20.5, 55.0, 3.0], [10.0, 10., 5.0], [0., 0., 0.], [100., 100., 50.]),
(:pseudovoigt, [20.5, 44.0, 3.0, 0.4], [5.0, 5., 5.0, 0.05], [0., 0., 0., 0.], [100., 100., 50., 1.])
]

result = fit_peaks2(x_fit, y_fit, peaks_info, sigma=noise, backend=:qGN, relax=6, maxiter=100)
```

You can plot the fit and components using:
```julia
result.plot_fit(components=true)
```

The fitted parameters can be printed using:
```julia
result.print_params(digits=4)
```

# References

Tarantola A., *Inverse Problem Theory and Methods for Model Parameter Estimation*, SIAM (2005), Chapter 3.

"""
function fit_peaks(x, y, peaks_info; sigma=nothing, backend=:Optim, relax=100, maxiter=100)
    param_indices = Dict()
    all_params = Float64[]
    all_prior_uncertainties = Float64[]
    all_lower_bounds = Float64[]
    all_upper_bounds = Float64[]
    current_idx = 1
    
    for (i, (peak_type, params, uncertainties, lower_bounds, upper_bounds)) in enumerate(peaks_info)
        n_params = peak_type == :pseudovoigt ? 4 : 3
        param_indices[i] = current_idx:(current_idx + n_params - 1)
        append!(all_params, params)
        append!(all_prior_uncertainties, uncertainties)
        append!(all_lower_bounds, lower_bounds)
        append!(all_upper_bounds, upper_bounds)
        current_idx += n_params
    end
    
    function model(params)
        # Use similar() instead of zeros() to match the type
        result = zero(x) * zero(params[1])  # This creates a vector of zeros with the right type
        
        for (i, (peak_type, _, _, _, _)) in enumerate(peaks_info)
            p = params[param_indices[i]]
            if peak_type == :gaussian
                result .+= gaussian(x, p)
            elseif peak_type == :lorentzian
                result .+= lorentzian(x, p)
            elseif peak_type == :pseudovoigt
                result .+= pseudovoigt(x, p)
            end
        end
        return result
    end

    
    function optimize_GN(params, prior_uncertainties; maxiter=maxiter, relax=relax)
      
        # data variance matrix
        CD = diagm(sigma.^2)
        ICD = inv(CD) # and its inverse

        # prior model variance matrix
        CM = diagm(prior_uncertainties.^2)
        ICM = inv(CM) # and its inverse
    
        #the model à priori, as a vector
        mprior = copy(params)
        mcurrent = copy(params) # initialize the algorithm at the à priori model values

        # loop
        iter = 1; # Just one variable to see the number of iteration

        while iter <= maxiter # So iterations will continue untill the difference of the chi2 is more than 1e-10. You can adjust that to your problem too

            ycalc = model(mcurrent) # model prediction
            G = ForwardDiff.jacobian(model, mcurrent) # derivatives
            GT = transpose(G)

            dataresiduals = ycalc-y # as the name variable says...
            modelresiduals = mcurrent-mprior # same

            # quasi Gauss-Newton algorithm
            # equation 3.51, page 69, Tarantola 2005
            gradient = GT*ICD*dataresiduals + CM\modelresiduals # The hessian gradient
            direction = (GT*ICD*G+ICM)\gradient # The direction
            mnew = mcurrent-direction./relax # with a relaxation parameter to avoid too large jumps, can be set to 1.

            mcurrent = deepcopy(mnew) # Here I set the current model to be equal to the new model
            iter = iter + 1
        end
        
        # derivatives at the optimum
        G = ForwardDiff.jacobian(model, mcurrent)
        GT = transpose(G)
        
        # formula 3.53, Tarantola 2005
        #CMpost = CM - CM*GT*inv(G*CM*GT + CD)*G*CM
        # FOR NOW IT DO NOT WORK... THE CMpost matrix is way too small.
        CMpost = inv(GT*ICD*G+ICM)
        errors = sqrt.(diag(CMpost))
            
        return mcurrent, errors, CMpost
    end
    
    # Define loss function with priors
    function loss(params)

        y_pred = model(params) # model predictions
        rss = sum((y .- y_pred).^2) # residual sum of squares
        
        # Add weighted residuals if weights are provided
        if sigma !== nothing
            rss = sum((y .- y_pred).^2/(sigma.^2))
        end

        # Add prior terms
        prior_loss = 0.0
        for (i, (_, init_params, uncertainties)) in enumerate(peaks_info)
            for j in 1:length(init_params)
                prior_loss += ((init_params[j] - params[param_indices[i][j]])^2)/(uncertainties[j]^2)
            end
        end
        
        return rss + prior_loss
    end
    
    # Perform the optimization
    result = optimize(loss, all_lower_bounds, all_upper_bounds, all_params, Fminbox(LBFGS()); autodiff = :forward)
    
    # Extract results
    if backend == :Optim
        fitted_params = Optim.minimizer(result)
        # Calculate errors (approximate)
        CMpost = nothing
        hessian_matrix = ForwardDiff.hessian(loss, fitted_params)
        errors = sqrt.(diag(inv(hessian_matrix)))
    elseif backend == :qGN
        fitted_params, errors, CMpost = optimize_GN(all_params, all_prior_uncertainties; maxiter=maxiter, relax=relax)
        hessian_matrix = ForwardDiff.hessian(loss, fitted_params)
        errors = sqrt.(diag(inv(hessian_matrix)))
    else
        error("Backend not implemented")
    end
    
    residuals = y - model(fitted_params)
    
    # Organize results by peak
    peak_results = []
    for (i, (peak_type, _, _)) in enumerate(peaks_info)
        push!(peak_results, (
            type = peak_type,
            params = fitted_params[param_indices[i]],
            errors = errors[param_indices[i]]
        ))
    end

    """
        plot_fit(; components=false)
    
        plot the data and the fit
    
    """
    function plot_fit(; components=false)
        p = scatter(x, y, label="Data", legend=:topleft)
        plot!(p, x, model(fitted_params), lw=2, label="Fit")
        
        if components
            for (i, peak) in enumerate(peak_results)
                component_params = fitted_params[param_indices[i]]
                if peak.type == :gaussian
                    plot!(p, x, gaussian(x, component_params), label="Gaussian $i")
                elseif peak.type == :lorentzian
                    plot!(p, x, lorentzian(x, component_params), label="Lorentzian $i")
                elseif peak.type == :pseudovoigt
                    plot!(p, x, pseudovoigt(x, component_params), label="PseudoVoigt $i")
                end
            end
        end
        
        return p
    end

    """
        print_params(; digits=4)
                
        print the parameters after the fit
    """
    function print_params(; digits=4)
        for (i, peak) in enumerate(peak_results)       
            println("Peak $i ($(peak.type)):")
            param_names = peak.type == :pseudovoigt ? 
                ["amplitude", "center", "width", "fraction"] : 
                ["amplitude", "center", "width"]
            
            for j in 1:length(peak.params)
                println("  $(param_names[j]): $(round(peak.params[j], digits=digits)) ± $(round(peak.errors[j], digits=digits))")
            end
        end
    end
                    
    # Return fit results and plotting function
    return (
        fit = result,
        peak_results = peak_results,
        residuals = residuals,
        hessian = hessian_matrix,
        plot_fit = plot_fit,
        print_params = print_params, 
        CMpost = CMpost
    )
end