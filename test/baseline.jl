using Spectra
using Test

# Dummy data
x = collect(0.:1.:100.)
y = collect(0.:1.:100.)
x2 = collect(-100.:0.1:100.)
y_poly4 = 1.0 .+ x2 .+ 0.25.*x2.^2 + 0.0001.*x2.^3 + 1.e-6.*x2.^4
y_peak = collect(0.:1.:100.) .+ 100.0.*exp.(-log(2).*((x.-50.0)./5.0).^2) .+ 0.2.*randn(size(x,1))
roi = [0.0 40;60. 100.0]

# testing baselines
y_fit_dierckx,bas_fit_dierckx = baseline(x[:],y,[0. 100.],"unispline",s=0.01)
@test isapprox(y_fit_dierckx,(y-y[:]),atol=1e-4)
@test isapprox(bas_fit_dierckx,y[:],atol=1e-4)

y_fit_linear,bas_fit_linear = baseline(x[:],y,[0. 100.],"poly",polynomial_order=2)
@test isapprox(y_fit_linear,(y-y[:]),atol=1e-4)
@test isapprox(bas_fit_linear,y[:],atol=1e-4)

y_fit_poly,bas_fit_poly = baseline(x2[:],y_poly4[:],[0. 100.],"poly",polynomial_order=4)
@test isapprox(y_fit_poly,(y_poly4[:]-y_poly4[:]),atol=1e-3)
@test isapprox(bas_fit_poly,y_poly4[:],atol=1e-3)

y_fit_als,bas_fit_als = baseline(x[:],y_peak,roi,"als",lam=10^5,p=0.001)
@test isapprox(y_fit_als,(y_peak-y[:]),atol=10.)
@test isapprox(bas_fit_als,y[:],atol=10.)

y_fit_arPLS,bas_fit_arPLS = baseline(x[:],y_peak,roi,"arPLS",lam=10^9,p=0.001)
@test isapprox(y_fit_arPLS,(y_peak-y[:]),atol=10.)
@test isapprox(bas_fit_arPLS,y[:],atol=10.)

#y_fit_whit,bas_fit_whit = baseline(x[:],y_peak[:],roi,"whittaker",lam=10^9)
#@test isapprox(y_fit_whit,(y_peak-y[:]),10.)
#@test isapprox(bas_fit_whit,y[:],10.)
