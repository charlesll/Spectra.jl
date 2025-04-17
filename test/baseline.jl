using Spectra
using Test

@testset "Baseline" begin
    # Dummy data
    x = collect(0.0:1.0:100.0)
    y = collect(0.0:1.0:100.0)
    x2 = collect(-100.0:0.1:100.0)
    y_poly4 = 1.0 .+ x2 .+ 0.25 .* x2 .^ 2 + 0.0001 .* x2 .^ 3 + 1.e-6 .* x2 .^ 4
    y_peak =
        y .+ 100.0 .* exp.(-log(2) .* ((x .- 50.0) ./ 5.0) .^ 2) .+ 0.2 .* randn(size(x, 1))
    roi = [0.0 40; 60.0 100.0]

    # Testing baselines
    # Polynomial baseline (order 2)
    y_fit, bas_ = baseline(x, y, roi=roi, method="polynomial", polynomial_order=2)
    @test isapprox(y, y_fit .+ bas_, atol=1e-4)

    # Polynomial baseline (order 4)
    y_fit, bas_ = baseline(x, y, roi=roi, method="polynomial", polynomial_order=4)
    @test isapprox(y, y_fit .+ bas_, atol=1e-4)

    # Dspline baseline
    y_fit, bas_ = baseline(x, y, roi=roi, method="Dspline", s=0.01)
    @test isapprox(y, y_fit .+ bas_, atol=1e-4)

    # gcvspline baseline
    y_fit, bas_ = baseline(x, y, roi=roi, method="gcvspline", s=1.0)
    @test isapprox(y, y_fit .+ bas_, atol=1e-4)

    # gcvspline baseline, no s provided
    y_fit, bas_ = baseline(x, y, roi=roi, method="gcvspline")
    @test isapprox(y, y_fit .+ bas_, atol=1e-4)

    # Whittaker smoother
    y_fit, bas_ = baseline(x, y_peak, roi=roi, method="whittaker", lambda=1.0e5)
    @test isapprox(y_peak, y_fit .+ bas_, atol=1e-4)

    # ALS baseline
    y_fit, bas_ = baseline(x, y_peak, method="als", lambda=1.0e5, p=0.001)
    @test isapprox(y_peak, y_fit .+ bas_, atol=1e-4)

    # arPLS baseline
    y_fit, bas_ = baseline(x, y_peak, method="arPLS", lambda=1.0e5, ratio=0.001)
    @test isapprox(y_peak, y_fit .+ bas_, atol=1e-4)

    # drPLS baseline
    y_fit, bas_ = baseline(x, y_peak, method="drPLS", lambda=1.0e5, ratio=0.1)
    @test isapprox(y_peak, y_fit .+ bas_, atol=1e-4)
end
