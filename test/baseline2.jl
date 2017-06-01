using Spectra
using Base.Test

# the x axis
x = collect(0:0.1:150)

# a scale factor that you can change
scale = 0.1

# a fake signal: perfect y
y_tot, y_peaks = gaussiennes([10;20.],[40.;60],[5.;15],x)
bkg, ~ = gaussiennes([10.],[50.],[50.],x)
y_perfect = scale.*(y_tot+bkg)

# we add noise: observed y
y = scale.*(y_tot+bkg + randn(size(y_tot,1)))

roi = [0 20.; 100. 150]

y1, bas1 = baseline(x,y,roi,"gcvspline",p=0.6)

rmse = sqrt(sum((bas1 - bkg).^2))./sum(bkg)

@test rmse <0.05 # error inferior to 5%, should actually be arround 0.026.