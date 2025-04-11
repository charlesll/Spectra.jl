using Spectra
using Test
using Random

@testset "peakmeas Tests" begin
    # Dummy data
    x = collect(0.:1.:100.)

    # Our peak is located at freq_th with a half-width of hwhm_th
    int_th = 10.
    freq_th = 40.
    hwhm_th = 10.

    # we now generate a perfect y as well as a noisy y
    y = int_th.*exp.(-log.(2) .*((x.-freq_th)./hwhm_th).^2)
    Random.seed!(42); # we fix the seed
    y_noise = y+randn((size(y,1),1))

    # we measured the width and frequency with the peakhw function
    int_mea1, freq_meas1, hwhm_meas1, centroid_mea1 = peakmeas(x,y,smoothing="no")
    #int_mea2, freq_meas2, hwhm_meas2, centroid_mea2 = peakmeas(x,vec(y_noise),method="savgol",window_length=5,polyorder=2) # with applying a Savitsky Golay filter

    # test of the perfect version
    @test isapprox(int_th,int_mea1,atol=1e-5)
    @test isapprox(freq_th,freq_meas1,atol=1e-5)
    @test isapprox(hwhm_th,hwhm_meas1,atol=1e-5)

    # test of the noisy version
    #@test isapprox(int_th,int_mea2,atol=1.)
    #@test isapprox(freq_th,freq_meas2,atol=1.)
    #@test isapprox(hwhm_th,hwhm_meas2,atol=1.)

    # Now we test 2 gaussians to test the centroid
    int2_th = 20.00000 # the two intensities cannot be fully equal, as this will raise an error.
    freq2_th = 60.
    hwhm2_th = 10.

    y2 = int2_th.*exp.(-log(2) .*((x.-freq2_th)./hwhm2_th).^2)
    y_tot = y + y2

    ~, ~, ~, centroid_mea3 = peakmeas(x,y_tot,smoothing="no")

    int_tot = int2_th + int_th

    centroid_th = int_th./int_tot.*freq_th .+ int2_th./int_tot.*freq2_th # this works because the peaks have the same hwhm
    @test isapprox(centroid_th,centroid_mea3,atol=1e-4)
end

@testset "centroid Tests" begin
    #### centroid test
    x = collect(0.:1.:100.)
    y_peak = gaussian(x, 1., 50., 10.)
    y_centroid = centroid(x, y_peak)
    @test isapprox(y_centroid, 50.0 ,atol=1e-4)

    y_centroid = centroid(x, y_peak, smoothing=true, method="gcvspline")
    @test isapprox(y_centroid, 50.0, atol=1e-4)
    
    y_centroids = centroid(x, [y_peak y_peak y_peak])
    @test isapprox(sum(y_centroids), 150.0, atol=1e-4)
    
    y_centroids = centroid([[x y_peak], [x y_peak], [x y_peak]])
    @test isapprox(sum(y_centroids), 150.0, atol=1e-4)
end

@testset "find_peaks Tests" begin
    # Test 1: Basic peak detection
    x = collect(1:0.1:100)
    y = gaussian(x, 1.0, 30.0, 3.0) + lorentzian(x, 1.0, 60.0, 6.0) + 0.01*randn(length(x))
    result = find_peaks(x, y; smoothing=false, window_size=15, min_height=0.5)
    
    @test length(result.peak_positions) > 0 #"Peak positions should be detected"
    @test length(result.peak_heights) == length(result.peak_positions) #"Peak heights should match detected positions"
    @test length(result.peak_prominences) == length(result.peak_positions) #"Peak prominences should match detected positions"
    @test length(result.peak_hwhms) == length(result.peak_positions) #"HWHMs should match detected positions"
    @test length(result.peak_centroids) == length(result.peak_positions) #"Centroids should match detected positions"
    @test isapprox(result.peak_positions, [30.0; 60.0], atol=1.0)

    # Test 2: Smoothing enabled
    result_smoothed = find_peaks(x, y; smoothing=true, window_size=10, min_height=0.2)
    @test length(result_smoothed.peak_positions) > 0 #"Peak positions should be detected with smoothing"

    # Test 3: Invalid input lengths
    x_invalid = collect(1:1.0:50)
    y_invalid = collect(1:1.0:100)
    
    @test_throws ErrorException find_peaks(x_invalid, y_invalid; smoothing=false, window_size=10)

    # Test 4: Invalid window size
    @test_throws ErrorException find_peaks(x, y; smoothing=false, window_size=-10)

    # Test 5: Plotting peaks (visual inspection)
    x = collect(1:0.1:100)
    y = gaussian(x, 1.0, 30.0, 3.0) + lorentzian(x, 1.0, 60.0, 6.0) + 0.01*randn(length(x))
    result_plot = find_peaks(x, y; smoothing=false, window_size=10)
    display(result_plot.plot_peaks)  # Ensure this produces a valid plot

end