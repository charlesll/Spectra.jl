using Spectra
using Base.Test

# Dummy data
x = collect(0.:1.:100.)

# Our peak is located at freq_th with a half-width of hwhm_th
freq_th = 50.
hwhm_th = 10.

# we now generate a perfect y as well as a noisy y
y = 10.0.*exp(-log(2) .*((x-freq_th)./hwhm_th).^2)
srand(42) # we fix the seed
y_noise = y+randn((size(y,1),1))

# we measured the width and frequency with the peakhw function
freq_meas1, hwhm_meas1 = peakhw(x,y)
freq_meas2, hwhm_meas2 = peakhw(x,y_noise,M=10,N=5) # with applying a Savitsky Golay filter


# test of the perfect version
@test_approx_eq_eps(freq_th,freq_meas1,1e-5)
@test_approx_eq_eps(hwhm_th,hwhm_meas1,1e-5)

# test of the noisy version
@test_approx_eq_eps(freq_th,freq_meas2,1.)
@test_approx_eq_eps(hwhm_th,hwhm_meas2,1.)
