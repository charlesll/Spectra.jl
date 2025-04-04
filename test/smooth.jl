using Spectra
using Test
using Statistics
using Random

@testset "Smooth" begin
    # the x axis
    x = collect(0:0.1:100)

    # a fake signal: perfect y
    y_tot, y_peaks = gaussiennes([10; 20.], [40.; 60], [5.; 15], x)
    y_perfect = vec(y_tot)

    # noise
    noise = randn(size(y_tot, 1))

    # observed y with added noise
    y = y_perfect .+ noise

    # Calculate the initial error (RMSE) between noisy signal and perfect signal
    initial_rmse = sqrt(mean((y_perfect - y).^2))

    # Loop over the methods and record the RMSE after smoothing
    methods_ = ["savgol", "whittaker", "gcvspline", "flat", "hanning", "hamming", "bartlett", "blackman"]
    for i in methods_
        y_smo = smooth(x, y, method=i)
        ese_method_ = sqrt(mean((y_perfect - y_smo).^2))  # RMSE after smoothing

        # Test that smoothing reduces the error compared to the initial noisy signal
        @test ese_method_ < initial_rmse
    end
end

@testset "Whittaker Smoother Tests" begin
    # Test case 1: Basic functionality with equally spaced x values
    @testset "Equally Spaced x" begin
        x = collect(1.0:0.1:10.0)
        y = sin.(x) .+ 0.1 .* randn(length(x))  # Noisy sine wave
        w = ones(length(x))                     # Equal weights
        lambda = 10.0                           # Smoothing parameter

        z = whittaker(x, y, w, lambda; d=2)

        # Check that output has the same length as input
        @test length(z) == length(y)

        # Check that smoothed values are closer to the true sine wave than the noisy input
        true_y = sin.(x)
        noisy_error = sqrt(mean((y .- true_y).^2))
        smoothed_error = sqrt(mean((z .- true_y).^2))
        @test smoothed_error < noisy_error
    end

    # Test case 2: Unequally spaced x values
    @testset "Unequally Spaced x" begin
        x = sort(rand(50) .* 10)                 # Randomly spaced x values in [0, 10]
        true_y = sin.(x)
        y = true_y .+ 0.1 .* randn(length(x))  # Noisy sine wave
        w = ones(length(x))                     # Equal weights
        lambda = 0.1                          # Smoothing parameter, ! its value depends on the length of x !!!
        z = whittaker(x, y, w, lambda; d=2)

        # Check that output has the same length as input
        @test length(z) == length(y)

        # Check that smoothed values are closer to the true sine wave than the noisy input
        noisy_error = sqrt(mean((y .- true_y).^2))
        smoothed_error = sqrt(mean((z .- true_y).^2))
        @test smoothed_error < noisy_error
    end

    # # Test case 3: Invalid inputs (error handling)
    # @testset "Error Handling" begin
    #     x = collect(1.0:10.0)
    #     y = sin.(x)
    #     w = ones(length(x))

    #     @test_throws ArgumentError whittaker(x[1:end-1], y, w, 10.0)   # Mismatched lengths of x and y
    #     @test_throws ArgumentError whittaker(x, y[1:end-1], w, 10.0)   # Mismatched lengths of y and w
    # end

    # Test case 4: Effect of smoothing parameter (lambda)
    @testset "Effect of Lambda" begin
        x = collect(1.0:0.1:10.0)
        y = sin.(x) .+ 0.5 .* randn(length(x))  # Noisy sine wave

        z_low_lambda = whittaker(x, y, ones(length(x)), 1.0; d=2)
        z_high_lambda = whittaker(x, y, ones(length(x)), 10000.0; d=2)

        # Low lambda should result in less smoothing (closer to noisy input)
        low_lambda_error = sqrt(mean((z_low_lambda .- y).^2))

        # High lambda should result in more smoothing (closer to a straight line or trend)
        high_lambda_error = sqrt(mean((z_high_lambda .- y).^2))

        @test low_lambda_error < high_lambda_error   # Low lambda retains more noise than high lambda
    end
end