using Spectra
using Test

# based on the example of the smooth function...

# the x axis
x = collect(0:0.1:100)

# a scale factor
scale = 0.01

# a fake signal: perfect y
y_tot, y_peaks = gaussiennes([10;20.],[40.;60],[5.;15],x)
y_perfect = scale.*y_tot

# we add noise: observed y
y = scale.*(y_tot + randn(size(y_tot,1)))

y_sv = smooth(x,y,filter=:SavitzkyGolay,M=15,N=2)
y_gcv = smooth(x,y,filter=:GCVSmoothedNSpline, ese_y = std(y))
y_whit = smooth(x,y,filter=:GCVSmoothedNSpline, lambda=10.0^1)

ese_noise = sum((y - y_perfect).^2)
ese_sg = sum((y_sv-y_perfect).^2)
ese_gcv = sum((y_gcv-y_perfect).^2)
ese_whit = sum((y_whit-y_perfect).^2)

@test ese_sg < ese_noise
@test ese_gcv < ese_noise
@test ese_whit < ese_noise
