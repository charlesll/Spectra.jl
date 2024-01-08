using Spectra
using Test

#### Dummy x
x = collect(0.1:0.1:100.)

#### gaussienne test
p = [10.,50.,5.]
y = p[1] .* exp.(-log(2) .* ((x.-p[2])./p[3]).^2)

y_calc, y_peaks = gaussiennes([p[1]],[p[2]],[p[3]],x)

@test isapprox(y,y_calc)

#### lorentzienne test
p = [10.,50.,5.]
y = p[1] ./ (1.0 .+ ((x.-p[2])./p[3]).^2)

y_calc, y_peaks = lorentziennes([p[1]],[p[2]],[p[3]],x)

@test isapprox(y,y_calc)

#### pearson7 test
a = [10.,50.,5., 5.]
y = a[1] ./ (1.0 .+ ((x.-a[2])./a[3]).^2 .* (2.0.^(1.0./a[4]) .- 1.0))

y_calc, y_peaks = pearson7([a[1]],[p[2]],[p[3]],[a[4]],x)

@test isapprox(y,y_calc)

#### pseudovoigt test
p = [10.,50.,5.,0.5]

y_lor = p[4] .* (p[1] ./ (1 .+ ((x.-p[2])./p[3]).^2))
y_gauss = (1.0 .-p[4]) .* (p[1] .* exp.(-log(2) .* ((x.-p[2])./p[3]).^2))
y = y_gauss + y_lor

y_calc, y_peaks = pseudovoigts([p[1]],[p[2]],[p[3]],[p[4]],x)
@test isapprox(y,y_calc)

#### normal distribution test
p = [1.,50.,5.]
y = p[1]./(p[3].*sqrt(2.0 .*pi)) .* exp.(-0.5.*((x.-p[2])./p[3]).^2)

y_calc, y_peaks = normal_dist([p[1]],[p[2]],[p[3]],x)
@test isapprox(y,y_calc)

#### centroid test
p = [10.,50.,5.]
y_calc, y_peaks = gaussiennes([p[1]],[p[2]],[p[3]],x)
y_centroid = centroid(x, y_calc)
println(size(y_centroid))
@test isapprox(50.,y_centroid[1,1])

#### normalisation test
y_norm1 = normalise(y_calc, method="intensity")
y_norm2 = normalise(y_calc, x=x, method="area")
y_norm3 = normalise(y_calc, method="minmax")

@test isapprox(1.0, maximum(y_norm1)[1,1])
@test isapprox(1.0, trapz(vec(x), vec(y_norm2))[1,1])
@test isapprox(1.0, maximum(y_norm3)[1,1])