using Spectra
using Test

# Dummy data
x = collect(0.:1.:100.)
y = collect(0.:1.:100.)
x2 = collect(-100.:0.1:100.)
y_poly4 = 1.0+ x2 + 0.25*x2.^2 + 0.0001.*x2.^3 + 1.e-6.*x2.^4
y_peak = collect(0.:1.:100.) + 100.0.*exp.(-log(2).*((x-50.0)/5.0).^2) + 0.2*randn(size(x,1))
roi = [0.0 40;60. 100.0]

# testing baselines
y_fit_gcvspl,bas_fit_gcvspl = baseline(x[:],y,[0. 100.],"gcvspline",p=0.01)
@test_approx_eq_eps(y_fit_gcvspl,(y-y[:]),1e-4)
@test_approx_eq_eps(bas_fit_gcvspl,y[:],1e-4)

y_fit_dierckx,bas_fit_dierckx = baseline(x[:],y,[0. 100.],"Dspline",p=0.01)
@test_approx_eq_eps(y_fit_dierckx,(y-y[:]),1e-4)
@test_approx_eq_eps(bas_fit_dierckx,y[:],1e-4)

y_fit_linear,bas_fit_linear = baseline(x[:],y,[0. 100.],"poly",p=2.0)
@test_approx_eq_eps(y_fit_linear,(y-y[:]),1e-4)
@test_approx_eq_eps(bas_fit_linear,y[:],1e-4)

y_fit_poly,bas_fit_poly = baseline(x2[:],y_poly4[:],[0. 100.],"poly",p=4.0)
@test_approx_eq_eps(y_fit_poly,(y_poly4[:]-y_poly4[:]),1e-3)
@test_approx_eq_eps(bas_fit_poly,y_poly4[:],1e-3)

#Tolerances for ML algo are high, it is just to test that they are running well...
y_fit_kr,bas_fit_kr = baseline(x[:],y[:],[0. 100.],"KRregression")
@test_approx_eq_eps(y_fit_kr,(y[:]-y[:]),1e-1)
@test_approx_eq_eps(bas_fit_kr,y[:],1e-1)

y_fit_svm,bas_fit_svm = baseline(x[:],y[:],[0. 100.],"SVMregression")
@test_approx_eq_eps(y_fit_svm,(y[:]-y[:]),10.)
@test_approx_eq_eps(bas_fit_svm,y[:],10.)

y_fit_als,bas_fit_als = baseline(x[:],y_peak,roi,"als",lambda=10^5,p=0.001)
@test_approx_eq_eps(y_fit_als,(y_peak-y[:]),10.)
@test_approx_eq_eps(bas_fit_als,y[:],10.)

y_fit_arPLS,bas_fit_arPLS = baseline(x[:],y_peak,roi,"arPLS",lambda=10^9,p=0.001)
@test_approx_eq_eps(y_fit_arPLS,(y_peak-y[:]),10.)
@test_approx_eq_eps(bas_fit_arPLS,y[:],10.)

y_fit_whit,bas_fit_whit = baseline(x[:],y_peak[:],roi,"whittaker",lambda=10^9)
@test_approx_eq_eps(y_fit_whit,(y_peak-y[:]),10.)
@test_approx_eq_eps(bas_fit_whit,y[:],10.)

