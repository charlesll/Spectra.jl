# Structure de contexte précalculé
struct FitContext
    x::Vector{Float64}
    peaks_info::Vector{Tuple}
    sigma::Vector{Float64}
    CM::Matrix{Float64}
    ICM::Matrix{Float64}
    CD::Matrix{Float64}
    ICD::Matrix{Float64}
    param_indices::Dict
    peak_functions::Dict
    mprior::Vector{Float64}
    mprior_sigma::Vector{Float64}
    all_lower_bounds::Vector{Float64}
    all_upper_bounds::Vector{Float64}
    model::Function  # Fonction prédéfinie
end

"""
    prepare_context(x, peaks_info, sigma=ones(length(x)))

Create a context for the fit.

`x` is the x-axis data
`peaks_info` is a vector of tuples with the peak type and parameters
    (peak_type, initial_params, uncertainties, lower_bounds, upper_bounds)
sigma is the noise in the data, default is ones(length(x))
"""
function prepare_context(x, peaks_info, sigma=ones(length(x)))
    # Calculs initiaux
    param_indices = Dict()
    mprior = Float64[]
    mprior_sigma = Float64[]
    all_lower_bounds = Float64[]
    all_upper_bounds = Float64[]
    current_idx = 1

    for (i, (peak_type, params, uncertainties, lower_bounds, upper_bounds)) in enumerate(peaks_info)
        n_params = peak_type ∈ (:pseudovoigt, :pearson7) ? 4 : 3
        param_indices[i] = current_idx:(current_idx + n_params - 1)
        append!(mprior, params)
        append!(mprior_sigma, uncertainties)
        append!(all_lower_bounds, lower_bounds)
        append!(all_upper_bounds, upper_bounds)
        current_idx += n_params
    end

    # Définition des matrices de covariance
    CM = diagm(mprior_sigma.^2)
    CD = diagm(sigma.^2)
    
    # peak dictionaries
    peak_functions = Dict(
        :gaussian => gaussian,
        :lorentzian => lorentzian,
        :pseudovoigt => pseudovoigt,
        :pearson7 => pearson7
    )

    # Définition permanente de la fonction model
    function model(params)
        result = zero(x) * zero(params[1])
        for (i, (peak_type, _, _, _, _)) in enumerate(peaks_info)
            p = params[param_indices[i]]
            result .+= peak_functions[peak_type](x, p)
        end
        return result
    end

    FitContext(
        x, peaks_info, sigma, CM, inv(CM), CD, inv(CD),
        param_indices, peak_functions, mprior,
        mprior_sigma, all_lower_bounds, all_upper_bounds,
        model
    )
end

"""
    fit_peaks(ctx::FitContext, y; backend=:qGN, relax=4, maxiter=100)

Perform a fit using the specified backend.

`ctx` is the context created by `prepare_context`
`y` is the data to fit
`backend` is the optimization backend to use, either :Optim or :qGN
`relax` is the relaxation factor for the step size (only used in :qGN)
`maxiter` is the maximum number of iterations (only used in :qGN)
"""
function fit_peaks(ctx::FitContext, y; backend=:Optim, relax=4, maxiter=100)
    
    # Perform the optimization
    if backend == :Optim
        fitted_params, CMpost, errors = fit_Optim(ctx, y)
    elseif backend == :qGN
        fitted_params, CMpost, errors = fit_qNewton(ctx, y, maxiter=maxiter, relax=relax)
    else
        error("Backend not implemented")
    end
    
    peak_results = get_peak_results(ctx, fitted_params, errors)
    # Retourner les résultats avec fermetures légères
    return (
        peak_results = peak_results,
        fitted_params = fitted_params,
        CMpost=CMpost,
        errors = errors,
        fit = ctx.model(fitted_params),
        residuals = y - ctx.model(fitted_params),
    )
end

"""
    fit_Optim(ctx::FitContext, y)
    
Perform a box constrained fit using the Optim package with the LBFGS algorithm.

The loss function combines a loss on data and on model prior (eq. 3.46 in Tarantola 2005)

`ctx` is the context created by `prepare_context`
`y` is the data to fit
"""
function fit_Optim(ctx::FitContext, y)
    # Définition des fermetures légères
    data_loss_dist = MvNormal(y, ctx.CD)
    mprior_loss_dist = MvNormal(ctx.mprior, ctx.CM)

    # loss by direct least-squares calculation
    # does not affect convergence, but do not use for Hessian calculation
    function ls_loss(params)
        data_loss = sum((ctx.model(params).-y).^2 ./ ctx.sigma.^2)
        mprior_loss = sum((ctx.mprior.-params).^2 ./ ctx.mprior_sigma.^2)
        out = data_loss + mprior_loss
        return out
    end
    function loss(params)
        data_loss = -logpdf(data_loss_dist,  ctx.model(params))
        mprior_loss = -logpdf(mprior_loss_dist, params)
        return data_loss + mprior_loss
    end

    # Perform the optimization
    result = optimize(
        loss,
        ctx.all_lower_bounds,
        ctx.all_upper_bounds,
        ctx.mprior,
        Fminbox(LBFGS());
        autodiff = :forward
    )
    fitted_params = Optim.minimizer(result)
    
    # Calculate errors (approximate)
    hessian_matrix = ForwardDiff.hessian(loss, fitted_params)
    CMpost = inv(hessian_matrix)
    errors = sqrt.(diag(CMpost))

    return fitted_params, CMpost, errors
end

"""
    fit_qNewton(ctx::FitContext, y; maxiter=100, relax=5)

Perform a quasi-Newton fit (Tarantola 2005, eq. 3.51)

`ctx` is the context created by `prepare_context`
`y` is the data to fit
`maxiter` is the maximum number of iterations
`relax` is the relaxation factor for the step size
"""
function fit_qNewton(ctx::FitContext, y; maxiter=100, relax=5)
    mcurrent = copy(ctx.mprior)
    
    for _ in 1:maxiter
        G = ForwardDiff.jacobian(ctx.model, mcurrent)
        GT = transpose(G)

        data_residuals = ctx.model(mcurrent) - y
        model_residuals = mcurrent - ctx.mprior

        gradient = GT * ctx.ICD * data_residuals + ctx.ICM * model_residuals
        direction = (GT * ctx.ICD * G + ctx.ICM) \ gradient
        mcurrent -= direction / relax
    end

    G = ForwardDiff.jacobian(ctx.model, mcurrent)
    GT = transpose(G)
    CMpost = ctx.CM - ctx.CM * GT * inv(G * ctx.CM * GT + ctx.CD) * G * ctx.CM
    return mcurrent, CMpost, sqrt.(diag(CMpost))
end

# Partie 4: Bootstrap optimisé
function bootstrap_fit(x, y, peaks_info, n=100; sigma=ones(length(x)))
    ctx = prepare_context(x, peaks_info, sigma)
    results = []
    
    Threads.@threads for _ in 1:n  # Parallelisation
        y_sample = y[rand(1:length(y), length(y))]
        push!(results, fit_peaks(ctx, y_sample))
    end
    
    results
end

"""
    plot_fit(ctx, params, y; components=false)

plot the data and the fit.

`ctx` is the context created by `prepare_context`
`params` is the vector of fitted parameters
`y` is the data to fit
`components` is a boolean to plot the components of the fit

"""
function plot_fit(ctx, params, y; components=false)
    p = scatter(ctx.x, y; label="Data", legend=:topleft)
    plot!(p, ctx.x, ctx.model(params); lw=2, label="Fit")

    if components
        for (i, (peak_type, _ , _ , _ , _ )) in enumerate(ctx.peaks_info)
            component_params = params[ctx.param_indices[i]]
            if peak_type == :gaussian
                plot!(p, ctx.x, gaussian(ctx.x, component_params); label="Gaussian $i")
            elseif peak_type == :lorentzian
                plot!(p, ctx.x, lorentzian(ctx.x, component_params); label="Lorentzian $i")
            elseif peak_type == :pseudovoigt
                plot!(p, ctx.x, pseudovoigt(ctx.x, component_params); label="PseudoVoigt $i")
            elseif peak_type == :pearson7
                plot!(p, ctx.x, pearson7(ctx.x, component_params); label="Pearson7 $i")
            end
        end
    end

    return p
end

"""
    print_params(peak_results; digits=4)
            
print the parameters after the fit, peak_results is a vector of Named Tuples
    with the results of the fit
"""
function print_params(peak_results; digits=2)
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
            # Arrondir paramètre et erreur
            rounded_param = round(peak.params[j], digits=digits)
            rounded_error_final = round(peak.errors[j], digits=digits)
            
            println(
                "  $(param_names[j]): $rounded_param ± $rounded_error_final",
            )
        
        end
    end
end

"""
    get_peak_results(ctx, fitted_params, fitted_params_errors)

Returns a vector of Named Tuples with the results of the fit

"""
function get_peak_results(ctx, fitted_params, fitted_params_errors)
    # Organize results by peak
    peak_results = []
    normalised_peak_areas = Float64[]
    for (i, (peak_type, _, _, _, _)) in enumerate(ctx.peaks_info)
        # get the parameters for the peak we are interested in
        params_ = fitted_params[ctx.param_indices[i]]
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
                peak_type = peak_type,
                params = params_,
                errors = fitted_params_errors[ctx.param_indices[i]],
                area = area_peak_out,
            ),
        )
    end
    return peak_results
end

"""
    bootstrap(x, y, sigma, peaks_info; nb_boot = 100, backend=:qGN, relax=5., maxiter=100)
    
Perform a bootstrap on the fit
"""
function bootstrap(x, y, sigma, peaks_info; nb_boot = 100, backend=:qGN, relax=5., maxiter=100)
    
    # preliminary fit
    ctx = prepare_context(x, peaks_info, sigma);
    result = fit_peaks(ctx, y, backend=backend, relax=relax, maxiter=maxiter)
    
    boot_params = zeros((length(result.fitted_params), nb_boot))
    for i in 1:nb_boot
        # get bootstrap sample
        x_boot, y_boot, ese_boot = bootsample(ctx.x, y, ese=ctx.sigma)
        ctx_boot = prepare_context(vec(x_boot), peaks_info, vec(ese_boot));
        # Perform the fit
        local result = fit_peaks(ctx_boot, vec(y_boot), backend=backend, relax=relax, maxiter=maxiter)
        boot_params[:,i] = result.fitted_params
    end
    return boot_params
end