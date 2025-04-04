using Spectra
using Test

@testset "Normalisation" begin
    # create data
    x = collect(0.:1.:100.)
    p = [10.,50.,5.,0.5]
    y_calc, y_peaks = pseudovoigts([p[1]],[p[2]],[p[3]],[p[4]],x)

    #### normalisation test
    y_norm1 = normalise(y_calc, method="intensity")
    y_norm2 = normalise(y_calc, x=x, method="area")
    y_norm3 = normalise(y_calc, method="minmax")

    y_norm_amm = normalise([y_norm1 y_norm2 y_norm3], method="intensity")

    @test isapprox(1.0, maximum(y_norm1)[1,1])
    @test isapprox(1.0, trapz(vec(x), vec(y_norm2))[1,1])
    @test isapprox(1.0, maximum(y_norm3)[1,1])
    @test isapprox(1.0, maximum(y_norm_amm)[1,1])
end


@testset "Resampling" begin
    x1 = [0., 2., 4.]
    y1 = 2.0*x1
    x_new = [1., 3.]

    # method 1
    y_new = resample(x1, y1, x_new)
    @test y_new == 2.0*x_new

    x2 = [1., 4., 8., 10.]
    y2 = [2.0*x2 x2]

    # method 2
    y_new = resample(x2, y2, x_new)
    @test y_new == [2.0*x_new x_new]

    # method 3
    y_new = resample([[x1 y1], [x1 y1]], x_new)
    @test y_new == [2.0*x_new 2.0*x_new]
end

@testset "extract_signal Tests" begin
    # Test case for single x-y input
    @testset "Single Signal" begin
        x = collect(1.0:10.0)
        y = x .^ 2
        roi = [3.0 5.0; 8.0 10.0]

        x_roi, y_roi, indices = extract_signal(x, y, roi)

        @test x_roi == [3.0, 4.0, 5.0, 8.0, 9.0, 10.0]
        @test y_roi == [9.0, 16.0, 25.0, 64.0, 81.0, 100.0]
        @test indices == [3, 4 ,5 ,8 ,9 ,10]
    end

    # Test case for spectrum matrix input
    @testset "Spectrum Matrix" begin
        x1 = collect(1.0:10.0)
        y1 = x1 .^ 2
        roi = [3.0 5.0; 8.0 10.0]

        x_roi, y_roi, indices = extract_signal([x1 y1], roi)

        @test x_roi == [3.0, 4.0, 5.0, 8.0, 9.0, 10.0]
        @test y_roi == [9.0, 16.0, 25.0, 64.0, 81.0, 100.0]
        @test indices == [3, 4 ,5 ,8 ,9 ,10]
    end

    # # Test case for multiple x-y inputs
    # @testset "Multiple Signals" begin
    #     x1 = collect(1.0:10.0)
    #     y1 = x1 .^ 2
    #     x2 = collect(11.0:20.0)
    #     y2 = x2 .^ 3
    #     roi = [3.0 5.0; 8.0 10.0]

    #     x_roi, y_roi, indices = extract_signal([x1, x2], [y1, y2], roi)

    #     @test x_roi == [[3.0, 4.0, 5.0], [8.0, 9.0, 10.0]]
    #     @test y_roi == [[9.0, 16.0, 25.0], [64.0, 81.0, 100.0]]
    #     @test indices == [[3,4 ,5], [8 ,9 ,10]]
    # end
end