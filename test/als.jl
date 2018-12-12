using Spectra
using Test

x = collect(50:1.0:500)

# gaussian peaks
p1 = 20.0 .* exp.(-log(2) .* ((x.-150.0)./15.0).^2)
p2 = 100.0 .* exp.(-log(2) .* ((x.-250.0)./5.0).^2)
p3 = 50.0 .* exp.(-log(2) .* ((x.-450.0)./1.0).^2)
p4 = 20.0 .* exp.(-log(2) .* ((x.-350.0)./30.0).^2)
p5 = 30.0 .* exp.(-log(2) .* ((x.-460.0)./5.0).^2)

# background: a large gaussian + linear 
bkg = 60.0 .* exp.(-log(2) .* ((x.-250.0)./200.0).^2) + 0.1.*x

#noise
noise = 2.0 * randn!(ones(size(x,1)))

#observation
y = p1 + p2 + p3 + p4 + p5 + noise +bkg

roi = [0. 100.] # provide whatever matrix

# for the ALS algorithm, 10^2-10^5 lambda and 0.001-0.1 p values are recommended
y_als, bas_als = baseline(x,y,roi,"als",p=0.04,lambda=10^6,niter=10)

rmse = sqrt.(sum((y - y_als).^2))./sum(y)

@test rmse <0.05 # error inferior to 5%, should actually be arround 0.026.