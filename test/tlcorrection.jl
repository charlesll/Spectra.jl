using Spectra
using Test

@testset "Basic Functionality" begin
    # Input data
    x = collect(100.0:10:1000)  # Raman shift in cm⁻¹
    y = rand(length(x))       # Random intensity values
    temperature_C = 25.0      # Temperature in Celsius
    laser_wavelength = 532.0  # Laser wavelength in nm

    # Run tlcorrection with default parameters
    x_out, y_corr, ese_corr = tlcorrection(x, y, temperature_C, laser_wavelength)

    # Check that output dimensions match input dimensions
    @test length(x_out) == length(x)
    @test length(y_corr) == length(y)
    @test length(ese_corr) == length(y)

    # Check that corrected values are not NaN or Inf
    @test all(.!isnan.(y_corr))
    @test all(.!isinf.(y_corr))
end

@testset "Long correction" begin
    # testing long correction function
    x_for_long = [20.0, 21.0, 22.0, 23.0, 24.0, 25.0]
    y_for_long = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    h_plank = 6.626070040e-34 #Plank constant
    k_bolt = 1.38064852e-23  #Boltzman constant
    c_light = 299792458.0 #speed of light in vaccum
    nu0 = 1.0 ./ (514.532) .* 1e9 #laser wavenumber at 514.532
    nu = 100.0 .* x_for_long # cm-1 to m-1
    T = 23.0+273.15 # the temperature in K

    x_long, long_res, eselong = tlcorrection([x_for_long y_for_long], 23.0, 514.532) # using the function
    t0 = nu0 .^ 3.0 .* nu ./ ((nu0 .- nu) .^ 4)
    t1 = 1 .- exp.(-h_plank .* c_light .* nu ./ (k_bolt .* T)) # c in m/s  : t1 dimensionless
    long_calc = y_for_long .* t0 .* t1 # pour les y
    long_calc = long_calc ./ trapz(x_for_long, long_calc) # area normalisation

    @test isapprox(long_res, long_calc)
    @test isapprox(x_for_long, x_long)
end

@testset "Correction Methods" begin
    x = collect(100.0:10:1000)
    y = rand(length(x))
    temperature_C = 25.0
    laser_wavelength = 532.0

    # Test "long" correction method
    x_long, y_long, ese_long = tlcorrection(
        x, y, temperature_C, laser_wavelength; correction="long"
    )
    @test all(.!isnan.(y_long))

    # Test "galeener" correction method
    x_galeener, y_galeener, ese_galeener = tlcorrection(
        x, y, temperature_C, laser_wavelength; correction="galeener"
    )
    @test all(.!isnan.(y_galeener))

    # Test "hehlen" correction method with density input
    x_hehlen, y_hehlen, ese_hehlen = tlcorrection(
        x, y, temperature_C, laser_wavelength; correction="hehlen", density=2200.0
    )
    @test all(.!isnan.(y_hehlen))
end

@testset "Normalisation Options" begin
    x = collect(100.0:10:1000)
    # a fake signal: perfect y
    y_peaks, y_tot = create_peaks(
        x,
        [
            Dict(:type => :gaussian, :amplitude => 10.0, :center => 40.0, :hwhm => 5.0),
            Dict(:type => :gaussian, :amplitude => 20.0, :center => 60.0, :hwhm => 15.0),
        ],
    )
    y = vec(y_tot)
    temperature_C = 25.0
    laser_wavelength = 532.0

    # Test area normalization
    _, y_area, _ = tlcorrection(x, y, temperature_C, laser_wavelength; normalisation="area")
    isapprox(sum(y_area), 1.0; atol=1e-6)

    # Test intensity normalization
    _, y_intensity, _ = tlcorrection(
        x, y, temperature_C, laser_wavelength; normalisation="intensity"
    )
    @test isapprox(maximum(y_intensity), 1.0; atol=1e-6)

    # Test minmax normalization
    _, y_minmax, _ = tlcorrection(
        x, y, temperature_C, laser_wavelength; normalisation="minmax"
    )
    @test isapprox(minimum(y_minmax), 0.0; atol=1e-6) &&
        isapprox(maximum(y_minmax), 1.0; atol=1e-6)

    # Test no normalization (output should match corrected values without scaling)
    _, y_no_norm, _ = tlcorrection(
        x, y, temperature_C, laser_wavelength; normalisation="no"
    )
end
