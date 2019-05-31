using Spectra
using Test
using Statistics
using Random

# based on the example of the smooth function...

# the x axis
x = collect(0:0.01:10)

# a fake signal: perfect y
y_perfect = x.^2

# fake noisy signal
y = y_perfect + randn(size(x,1))*2

y_sv = smooth(x,y,method="savgol",window_length=9,polyorder=3)
y_whit = smooth(x,y,method="whittaker", Lambda=10.0^1)

ese_noise = sum((y .- y_perfect).^2)
ese_sg = sum((y_sv.-y_perfect).^2)
ese_whit = sum((y_whit.-y_perfect).^2)

@test ese_sg < ese_noise
@test ese_whit < ese_noise
