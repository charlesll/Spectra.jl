using Spectra
using Test

@testset "Peak shapes and creation Tests" begin
    x = collect(1.0:0.5:100.0)
    a = [10.0, 50.0, 5.0, 0.5]

    gauss_ = a[1] .* exp.(-log(2) .* ((x .- a[2]) ./ a[3]) .^ 2)
    lor_ = a[1] ./ (1.0 .+ ((x .- a[2]) ./ a[3]) .^ 2)
    pv_ = a[4]*lor_ + (1-a[4])*gauss_
    pears_ =
        y =
            a[1] ./
            (1.0 .+ ((x .- a[2]) ./ a[3]) .^ 2 .* (2.0 .^ (1.0 ./ a[4]) .- 1.0)) .^ a[4]

    @test isapprox(gauss_, gaussian(x, a[1:3]))
    @test isapprox(lor_, lorentzian(x, a[1:3]))
    @test isapprox(pv_, pseudovoigt(x, a))
    @test isapprox(pears_, pearson7(x, a))
    @test isapprox(gauss_, gaussian(x, a[1:3]...))
    @test isapprox(lor_, lorentzian(x, a[1:3]...))
    @test isapprox(pv_, pseudovoigt(x, a...))
    @test isapprox(pears_, pearson7(x, a...))

    peak_infos = [
        Dict(:type => :gaussian, :amplitude => a[1], :center => a[2], :hwhm => a[3]),
        Dict(:type => :lorentzian, :amplitude => a[1], :center => a[2], :hwhm => a[3]),
        Dict(
            :type => :pseudovoigt,
            :amplitude => a[1],
            :center => a[2],
            :hwhm => a[3],
            :lorentzian_fraction => a[4],
        ),
        Dict(
            :type => :pearson7,
            :amplitude => a[1],
            :center => a[2],
            :hwhm => a[3],
            :exponent => a[4],
        ),
    ]
    peak_calculated = [gauss_, lor_, pv_, pears_]
    peaks, total = create_peaks(x, peak_infos)
    for i in 1:size(peaks, 2)
        @test isapprox(peak_calculated[i], peaks[:, i])
    end

    #### normal distribution test
    p = [1.0, 50.0, 5.0]
    y = p[1] ./ (p[3] .* sqrt(2.0 .* pi)) .* exp.(-0.5 .* ((x .- p[2]) ./ p[3]) .^ 2)

    @test isapprox(y, normal_dist(x, p...))
end
