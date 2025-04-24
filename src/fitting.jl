"""
FitContext type containing all fitting context
-`x::Vector{Float64}`: x data vector.
-`y::Vector{Float64}`: y data vector.
-`peaks_info::Vector{Tuple}`: contains the informations of the peaks of the model.
-`sigma::Vector{Float64}`: error on y.
-`CD::Matrix{Float64}`: data covariance matrix.
-`ICD::Matrix{Float64}`: inverse of the data covariance matrix.
-`mprior::Vector{Float64}`: prior model parameters.
-`mprior_sigma::Vector{Float64}`: uncertainties on prior model parameters.
-`CM::Matrix{Float64}`: model covariance matrix (prior).
-`ICM::Matrix{Float64}`: inverse of the model covariance matrix (prior).
-`all_lower_bounds::Vector{Float64}`: lower boundaries for parameters.
-`all_upper_bounds::Vector{Float64}`: upper boundaries for parameters.
-`model::Function`
"""
struct FitContext
    x::Vector{Float64}
    y::Vector{Float64}
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
FitResult type containing all fitting results and metadata
`context::FitContext`: Fit Context.
`peak_results::Vector`: Vector of Named Tuples containing parameters and areas for each peak, with errors.
`params::Vector{Float64}`: Vector of parameters.
`covariance::Matrix{Float64}`: Covariance matrix at the optimal point ``C_{M, post}``.
`errors::Vector{Float64}`: Uncertainties on parameters calculated as ``sqrt.(diag(covariance))``.
`y_calc::Vector{Float64}`: Calculated y values, using `params`.
`residuals::Vector{Float64}`: Difference between calculated and input y values.
"""
struct FitResult
    context::FitContext
    peak_results::Vector
    params::Vector{Float64}
    covariance::Matrix{Float64}
    errors::Vector{Float64}
    y_calc::Vector{Float64}
    residuals::Vector{Float64}
    #timestamp::DateTime
    #convergence::Symbol
end
"""
    prepare_context(x, y, peaks_info, sigma=ones(length(x)))

Create a precomputed context for peak fitting operations.

# Arguments
- `x::Vector{Float64}`: Independent variable values
- `y::Vector{Float64}`: Observations
- `peaks_info::Vector{Tuple}`: Peak specifications containing:
  - `Symbol`: Peak type (`:gaussian`, `:lorentzian`, `:pseudovoigt`, `:pearson7`)
  - `Vector{Float64}`: Initial parameters
  - `Vector{Float64}`: Parameter uncertainties
  - `Vector{Float64}`: Lower bounds
  - `Vector{Float64}`: Upper bounds
- `sigma::Vector{Float64}`: Data uncertainties (default: ones)

# Returns
- `FitContext`: Precomputed structure containing matrices, indices, and model function

# Examples
'''julia
peaks_info = [
(:gaussian, [1.0, 0.0, 0.5], [0.1, 0.05, 0.1], [-Inf, -1.0, 0.1], [Inf, 1.0, 2.0]),
(:lorentzian, [0.5, 0.2, 0.3], [0.1, 0.05, 0.1], [0.0, -0.5, 0.1], [2.0, 0.5, 1.0])
]
ctx = prepare_context(x, peaks_info, sigma)
'''
"""
function prepare_context(x, y, peaks_info, sigma=ones(length(x)))
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
        x, y, peaks_info, sigma, CM, inv(CM), CD, inv(CD),
        param_indices, peak_functions, mprior,
        mprior_sigma, all_lower_bounds, all_upper_bounds,
        model
    )
end

"""
    fit_peaks(ctx::FitContext; backend=:qGN, relax=4, maxiter=100)

Perform peak fitting using specified optimization backend.

# Arguments
- `ctx::FitContext`: Precomputed context from `prepare_context`
- `backend::Symbol`: Optimization method (`:Optim` or `:qGN`)
- `relax::Real`: Step relaxation factor (qGN only)
- `maxiter::Int`: Maximum iterations (qGN only)

# Returns
- `FitResult`: Structured results containing:
  - `context::FitContext`: Fit context
  - `peak_results::Vector`: Peak results with uncertainties
  - `params::Vector{Float64}`: Peak parameters
  - `covariance::Matrix{Float64}`: Covariance matrix
  - `errors::Vector{Float64}`: 1-sigma standard errors on parameters
  - `y_calc::Vector{Float64}`: Model predictions
  - `residuals::Vector{Float64}`: Residuals

# Examples
```julia
x = 1:1.0:100
y = gaussian(x, 1., 50.0, 5.0) .+ 0.1 * randn(length(x))
peaks_info = [
    (:gaussian, [1.0, 50.0, 5.0], [0.1, 0.05, 0.1], [0.0, 40.0, 0.1], [2.0, 60.0, 10.0])
]
ctx = prepare_context(x, y, peaks_info)
result = fit_peaks(ctx; backend=:qGN, relax=4, maxiter=100)
plot_fit(ctx; result=result.peak_results)
```
"""
function fit_peaks(ctx::FitContext; backend=:Optim, relax=4, maxiter=100)
    
    # Perform the optimization
    if backend == :Optim
        fitted_params, CMpost, errors = fit_Optim(ctx)
    elseif backend == :qGN
        fitted_params, CMpost, errors = fit_qNewton(ctx, maxiter=maxiter, relax=relax)
    else
        error("Backend not implemented")
    end
    
    peak_results = get_peak_results(ctx, fitted_params, CMpost)
    # Retourner les résultats avec fermetures légères
    out = FitResult(ctx, 
        peak_results, 
        fitted_params, 
        CMpost, errors, 
        ctx.model(fitted_params), 
        ctx.y - ctx.model(fitted_params))
    return out
end

"""
    fit_Optim(ctx::FitContext)
    
Perform a box constrained fit using the Optim package with the LBFGS algorithm.

The loss function combines a loss on data and on model prior (eq. 3.46 in Tarantola 2005)

# Arguments
-`ctx`: context created by `prepare_context`

# Returns
- `fitted_params`: fitted parameters
- `CMpost`: covariance matrix of the fitted parameters
- `sqrt.(diag(CMpost))`: errors of the fitted parameters
"""
function fit_Optim(ctx::FitContext)
    # Définition des fermetures légères
    data_loss_dist = MvNormal(ctx.y, ctx.CD)
    mprior_loss_dist = MvNormal(ctx.mprior, ctx.CM)

    # loss by direct least-squares calculation
    # does not affect convergence, but do not use for Hessian calculation
    function ls_loss(params)
        data_loss = sum((ctx.model(params).-ctx.y).^2 ./ ctx.sigma.^2)
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
    dfc = TwiceDifferentiableConstraints(ctx.all_lower_bounds, ctx.all_upper_bounds)
    result = optimize(
        ls_loss, # we use the least-squares loss for the optimization
        dfc,
        ctx.mprior,
        IPNewton();
        autodiff = :forward
    )
    fitted_params = Optim.minimizer(result)
    
    # Calculate errors, this time we use the full (scaled) loss = log likelihood + log prior
    hessian_matrix = ForwardDiff.hessian(loss, fitted_params)
    CMpost = inv(hessian_matrix)
    errors = sqrt.(diag(CMpost))

    return fitted_params, CMpost, errors
end

"""
    fit_qNewton(ctx::FitContext; maxiter=100, relax=5)

Perform a quasi-Newton fit (Tarantola 2005, eq. 3.51)

# Arguments
-`ctx`: context created by `prepare_context`
-`maxiter`: maximum number of iterations
-`relax`: relaxation factor for the step size

# Returns
- `mcurrent`: fitted parameters
- `CMpost`: covariance matrix of the fitted parameters
- `sqrt.(diag(CMpost))`: errors of the fitted parameters
"""
function fit_qNewton(ctx::FitContext; maxiter=100, relax=5)
    mcurrent = copy(ctx.mprior)
    
    for _ in 1:maxiter
        G = ForwardDiff.jacobian(ctx.model, mcurrent)
        GT = transpose(G)

        data_residuals = ctx.model(mcurrent) - ctx.y
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

"""
    plot_fit(ctx, peak_results=nothing; xlabel="X", ylabel="Y", title="Model adjustement")

return a plot of the data and the fit given a `ctx` context created by `prepare_context`, 
and a `result` generated by `fit_peaks` or `get_fit_results`.
If result is not provided, the prior is plotted.

"""
function plot_fit(ctx, peak_results=nothing; xlabel="X", ylabel="Y", title="Model adjustement")
    p1 = scatter(ctx.x, ctx.y; ms=1, mc=:black,
    label="Data", 
    ylabel=ylabel,
    title=title,
    legend=:topleft)

    if peak_results == nothing
        println("No result provided, plotting the prior...")
        iter_NamedVector = ctx.peaks_info
        plot!(p1, ctx.x, ctx.model(ctx.mprior); lw=1, lc=:red, label="Prior model")
        p2 = plot(ctx.x, ctx.model(ctx.mprior)-ctx.y, xlabel=xlabel, ylabel="Residuals", legend=false)
    else
        iter_NamedVector = peak_results
        vectorized_params = []
        for (_, params, _) in iter_NamedVector
            append!(vectorized_params, params)
        end
        y_calc = Measurements.value.(ctx.model(vectorized_params))
        plot!(p1, ctx.x, y_calc; lw=1, lc=:red, label="Adjusted model")
        p2 = plot(ctx.x, y_calc-ctx.y, xlabel=xlabel, ylabel="Residuals", legend=false)
    end

    for (i, (peak_type, params , _ )) in enumerate(iter_NamedVector)
        if peak_type == :gaussian
            plot!(p1, ctx.x, Measurements.value.(gaussian(ctx.x, params)); label="Gaussian $i")
        elseif peak_type == :lorentzian
            plot!(p1, ctx.x, Measurements.value.(lorentzian(ctx.x, params)); label="Lorentzian $i")
        elseif peak_type == :pseudovoigt
            plot!(p1, ctx.x, Measurements.value.(pseudovoigt(ctx.x, params)); label="PseudoVoigt $i")
        elseif peak_type == :pearson7
            plot!(p1, ctx.x, Measurements.value.(pearson7(ctx.x, params)); label="Pearson7 $i")
        end
    end

    p = plot(p1, p2, layout=grid(2, 1, heights=[0.8, 0.2]), size=(800, 600))
    return p
end

"""
    print_params(peak_results)
            
print the parameters after the fit, peak_results is a vector of Named Tuples 
    (generated by get_peak_results) with the results of the fit.
"""
function print_params(peak_results)
    for (i, peak) in enumerate(peak_results)
        println("Peak $i ($(peak.peak_type)):")
        if peak.peak_type == :pseudovoigt
            param_names = ["amplitude", "center", "width", "fraction"]
        elseif peak.peak_type == :pearson7
            param_names = ["amplitude", "center", "width", "exponent"]
        else
            param_names = ["amplitude", "center", "width"]
        end
        # printing
        for j in 1:length(peak.params)            
            println("  $(param_names[j]): $(peak.params[j])")
        end
        # printing area
        println("  area: $(peak.area)")
    end
end

"""
    get_peak_results(ctx, params, CMpost)

Returns a vector of Named Tuples with the results of the fit. 
Values are tied to their errors via the covariance matrix CMpost.

"""
function get_peak_results(ctx, params, CMpost)
    
    # first we create a new params array containing the errors
    params_with_ese = Measurements.correlated_values(params, CMpost)
    
    # then we organize results by peak
    peak_results = []
    for (i, (peak_type, _, _, _, _)) in enumerate(ctx.peaks_info)
        # get the parameters for the peak we are interested in
        params_ = params_with_ese[ctx.param_indices[i]]
        amplitude_ = params_[1]
        hwhm_ = params_[3]
        lorentzian_fraction_ = peak_type == :pseudovoigt ? params_[4] : nothing
        exponent_ = peak_type == :pearson7 ? params_[4] : nothing

        # calculate the area and store it in a vector of the normalised peak areas
        area_ = area_peaks(
            peak_type,
            amplitude_,
            hwhm_;
            lorentzian_fraction=lorentzian_fraction_,
            exponent=exponent_,
        )

        # store the peak results also in peak_results Named Tuple
        push!(
            peak_results,
            (
                peak_type = peak_type,
                params = params_,
                area = area_,
            ),
        )
    end
    return peak_results
end

"""
    bootstrap(x, y, sigma, peaks_info; nb_boot = 100, backend=:qGN, relax=5., maxiter=100)
    
Perform a bootstrap on the fit.

# Arguments
-`x`: x-axis data
-`y`: y-axis data
-`sigma`: data noise
-`peaks_info`: vector of tuples with the peak type and parameters
    (peak_type, initial_params, uncertainties, lower_bounds, upper_bounds)
-`nb_boot`: number of bootstrap samples
-`backend`: optimization backend, either :Optim or :qGN
-`relax`: relaxation factor for the step size (only used in :qGN)
-`maxiter`: maximum number of iterations (only used in :qGN)

# Returns
- a matrix of size (number of parameters, number of bootstrap samples) with the fitted parameters
- a vector of Named Tuples with the results of the fit. Values are tied to their errors via the covariance matrix CMpost.
    (peak_type, initial_params, uncertainties, lower_bounds, upper_bounds)  
"""
function bootstrap(x, y, sigma, peaks_info; nb_boot = 100, backend=:qGN, relax=5., maxiter=100)
    
    # preliminary fit
    ctx = prepare_context(x, y, peaks_info, sigma);
    result = fit_peaks(ctx, backend=backend, relax=relax, maxiter=maxiter)
    
    boot_params = zeros((length(result.params), nb_boot))
    for i in 1:nb_boot
        # get bootstrap sample
        x_boot, y_boot, ese_boot = bootsample(ctx.x, ctx.y, ese=ctx.sigma)
        ctx_boot = prepare_context(vec(x_boot), vec(y_boot), peaks_info, vec(ese_boot));
        # Perform the fit
        local result = fit_peaks(ctx_boot, backend=backend, relax=relax, maxiter=maxiter)
        boot_params[:,i] = result.params
    end

    # We get the results by passing the median and standard deviation of the bootstrap samples.
    # Remember: the samples are along the Array columns.
    boot_results = get_peak_results(ctx, vec(Statistics.median(boot_params, dims=2)), Statistics.cov(boot_params, dims=2))

    return boot_params, boot_results
end