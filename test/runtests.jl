using Spectra
using Base.Test

# Dummy data
x = collect(0.:1.:100.)
y = collect(0.:1.:100.)
x2 = collect(-100.:0.1:100.)
y_poly4 = 1.0+ x2 + 0.25*x2.^2 + 0.0001.*x2.^3 + 1.e-6.*x2.^4
roi = [0.0 100.0]

# testing baselines
y_fit_gcvspl,bas_fit_gcvspl = baseline(x[:],y[:],roi,"gcvspline",p=1.0)
@test_approx_eq_eps(y_fit_gcvspl,(y[:]-y[:]),1e-10)
@test_approx_eq_eps(bas_fit_gcvspl,y[:],1e-10)

y_fit_dierckx,bas_fit_dierckx = baseline(x[:],y[:],roi,"Dspline",p=1.0)
@test_approx_eq_eps(y_fit_dierckx,(y[:]-y[:]),1e-10)
@test_approx_eq_eps(bas_fit_dierckx,y[:],1e-10)

y_fit_linear,bas_fit_linear = baseline(x[:],y[:],roi,"poly",p=2.0)
@test_approx_eq_eps(y_fit_linear,(y[:]-y[:]),1e-10)
@test_approx_eq_eps(bas_fit_linear,y[:],1e-10)

y_fit_poly,bas_fit_poly = baseline(x2[:],y_poly4[:],roi,"poly",p=4.0)
@test_approx_eq_eps(y_fit_poly,(y_poly4[:]-y_poly4[:]),1e-3)
@test_approx_eq_eps(bas_fit_poly,y_poly4[:],1e-3)

#Tolerances for ML algo are high, it is just to test that they are running well...
y_fit_kr,bas_fit_kr = baseline(x[:],y[:],roi,"KRregression")
@test_approx_eq_eps(y_fit_kr,(y[:]-y[:]),1e-1)
@test_approx_eq_eps(bas_fit_kr,y[:],1e-1)

y_fit_svm,bas_fit_svm = baseline(x[:],y[:],roi,"SVMregression")
@test_approx_eq_eps(y_fit_svm,(y[:]-y[:]),10.)
@test_approx_eq_eps(bas_fit_svm,y[:],10.)

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

#figure()
#plot(x2,y_poly4,"o")
#plot(x2,bas_fit_poly,"-")
