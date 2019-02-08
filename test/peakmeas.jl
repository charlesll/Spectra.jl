using Spectra
using Test
using Random

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
int_mea2, freq_meas2, hwhm_meas2, centroid_mea2 = peakmeas(x,y_noise,M=9,N=5) # with applying a Savitsky Golay filter

# test of the perfect version
@test isapprox(int_th,int_mea1,1e-5)
@test isapprox(freq_th,freq_meas1,1e-5)
@test isapprox(hwhm_th,hwhm_meas1,1e-5)

# test of the noisy version
@test isapprox(int_th,int_mea2,1.)
@test isapprox(freq_th,freq_meas2,1.)
@test isapprox(hwhm_th,hwhm_meas2,1.)

# Now we test 2 gaussians to test the centroid
int2_th = 20.00000 # the two intensities cannot be fully equal, as this will raise an error.
freq2_th = 60.
hwhm2_th = 10.

y2 = int2_th.*exp.(-log(2) .*((x-freq2_th)./hwhm2_th).^2)
y_tot = y + y2

~, ~, ~, centroid_mea3 = peakmeas(x,y_tot,smoothing="no")

int_tot = int2_th + int_th

centroid_th = int_th./int_tot.*freq_th + int2_th./int_tot.*freq2_th # this works because the peaks have the same hwhm
@test isapprox(centroid_th,centroid_mea3,1e-4)
