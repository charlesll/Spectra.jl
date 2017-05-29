using Spectra
using Base.Test

# testing long correction function
x_for_long = [20.,21.,22.,23.,24.,25.]
y_for_long = [1.0,1.0,1.0,1.0,1.0,1.0]
h_plank = 6.626070040e-34 #Plank constant
k_bolt = 1.38064852e-23  #Boltzman constant
c_light = 299792458. #speed of light in vaccum
nu0 = 1.0./(514.532).*1e9 #laser wavenumber at 514.532
nu = 100.0.*x_for_long # cm-1 to m-1
T = 23.0+273.15 # the temperature in K

x_long,long_res,eselong = tlcorrection([x_for_long y_for_long],23.0,514.532) # using the function
t0 = nu0.^3.*nu./((nu0-nu).^4)
t1= 1 - exp(-h_plank.*c_light.*nu./(k_bolt.*T)) # c in m/s  : t1 dimensionless
long_calc= y_for_long.*t0.*t1 # pour les y
long_calc = long_calc./trapz(x_for_long,long_calc) # area normalisation

@test_approx_eq(long_res,long_calc)
@test_approx_eq(x_for_long,x_long)
