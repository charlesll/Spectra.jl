using Spectra
using Test

# Dummy data
x = collect(-100.:1.:100.)

p = [10.,1.0,0.25]
y = p[1]+ p[2].*x + p[3].*x.^2

# using the poly and polyfit functions
y_calc = poly(p,x)
coef_calc = polyfit(x,y,2)

@test_approx_eq_eps(y,y_calc,1e-5)
@test_approx_eq_eps(p,coef_calc,1e-5)

