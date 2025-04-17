using Test, Spectra

@testset "Optim" begin
    # Define peak information with initial parameters and uncertainties
    x_fit = sort(rand(1000)*100)
    noise = randn(length(x_fit))*0.3

    # Gaussian case
    peak = gaussian(x_fit, [10.5, 30.0, 3.0])
    peaks_info = [
    # (type, initial_params, uncertainties, lower_bounds, upper_bounds)
        (
        :gaussian,
        [10.5, 30.0, 3.0],
        [100.0, 100.0, 100.0],
        [0.0, 0.0, 0.0],
        [Inf, Inf, Inf],
    ),]

    ctx = prepare_context(x_fit, peaks_info, noise)
    result = fit_peaks(ctx, peak + noise, backend=:Optim)
    @test isapprox(result.fit, peak, rtol=0.1)

    # Lorentzian
    peak = lorentzian(x_fit, [10.5, 30.0, 3.0])
    peaks_info = [
    # (type, initial_params, uncertainties, lower_bounds, upper_bounds)
        (
        :lorentzian,
        [10.5, 30.0, 3.0],
        [100.0, 100.0, 100.0],
        [0.0, 0.0, 0.0],
        [Inf, Inf, Inf],
    ),
    ]

    ctx = prepare_context(x_fit, peaks_info, noise)
    result = fit_peaks(ctx, peak + noise, backend=:Optim)
    @test isapprox(result.fit, peak, rtol=0.1)

    # pseudovoigt
    peak = pseudovoigt(x_fit, [10.0, 35.0, 10.0, 0.5])
    peaks_info = [
    # (type, initial_params, uncertainties, lower_bounds, upper_bounds)
        (
        :pseudovoigt,
        [10.5, 30.0, 3.0, 0.5],
        [100.0, 100.0, 100.0, 1.0],
        [0.0, 0.0, 0.0, 0.0],
        [Inf, Inf, Inf, 1.0],
    ),
    ]

    ctx = prepare_context(x_fit, peaks_info, noise)
    result = fit_peaks(ctx, peak + noise, backend=:Optim)
    @test isapprox(result.fit, peak, rtol=0.1)
end

@testset "Optim, multiple peaks" begin
    x_fit = sort(rand(1000)*100)
    y_fit_perfect = (
        gaussian(x_fit, [10.0, 35.0, 10.0]) +
        lorentzian(x_fit, [15.0, 55.0, 3.0]) +
        pseudovoigt(x_fit, [20.0, 45.0, 2.0, 0.4])
    )
    noise = randn(length(x_fit))*0.3
    y_fit = y_fit_perfect + noise

    peaks_info = [
        # (type, initial_params, uncertainties, lower_bounds, upper_bounds)
        (:gaussian, [10.5, 30.0, 3.0], [5.0, 10.0, 5.0], [0.0, 0.0, 0.0], [Inf, Inf, 50.0]),
        (
            :lorentzian,
            [15.5, 55.0, 3.1],
            [5.0, 10.0, 5.0],
            [0.0, 0.0, 0.0],
            [100.0, 100.0, 50.0],
        ),
        (
            :pseudovoigt,
            [20.5, 44.0, 3.0, 0.4],
            [5.0, 10.0, 5.0, 0.02],
            [0.0, 0.0, 0.0, 0.0],
            [100.0, 100.0, 50.0, 1.0],
        ),
    ]

    ctx = prepare_context(x_fit, peaks_info, noise)
    result = fit_peaks(ctx, y_fit, backend=:Optim)

    @test isapprox(result.fit, y_fit_perfect, rtol=0.1)
end

@testset "quasi Gauss Newton, multiple peaks" begin
    x_fit = sort(rand(1000)*100)
    y_fit_perfect = (
        gaussian(x_fit, [10.0, 35.0, 10.0]) +
        lorentzian(x_fit, [15.0, 55.0, 3.0]) +
        pseudovoigt(x_fit, [20.0, 45.0, 2.0, 0.4])
    )
    noise = randn(length(x_fit))*0.3
    y_fit = y_fit_perfect + noise

    peaks_info = [
        # (type, initial_params, uncertainties, lower_bounds, upper_bounds)
        (:gaussian, [10.5, 30.0, 3.0], [5.0, 10.0, 5.0], [0.0, 0.0, 0.0], [Inf, Inf, 50.0]),
        (
            :lorentzian,
            [20.5, 55.0, 3.0],
            [10.0, 10.0, 5.0],
            [0.0, 0.0, 0.0],
            [100.0, 100.0, 50.0],
        ),
        (
            :pseudovoigt,
            [20.5, 44.0, 3.0, 0.4],
            [5.0, 5.0, 5.0, 0.05],
            [0.0, 0.0, 0.0, 0.0],
            [100.0, 100.0, 50.0, 1.0],
        ),
    ]

    ctx = prepare_context(x_fit, peaks_info, noise)
    result = fit_peaks(ctx, y_fit, backend=:qGN, relax=100, maxiter=1000)
    
    @test isapprox(result.fit, y_fit_perfect, rtol=0.1)
end
