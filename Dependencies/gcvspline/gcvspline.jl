using PyCall
unshift!(PyVector(pyimport("sys")["path"]), "")

@pyimport gcvspline

function Gspline(xfit::Array{Float64},yfit::Array{Float64},esefit::Array{Float64},xtarget::Array{Float64},smoothing::Float64;VAL = esefit.^2,splor=2,splm=3,nc = size(yfit,1))
    c, wk, ier = gcvspline.gcvspline(xfit,yfit,smoothing.*ese,VAL,splmode=splm, splorder=splor,NC=nc) # gcvspl with mode 3 and smooth factor
    out = gcvspline.splderivative(x,xfit,c,splorder=splor, L = 1,IDER = 0) 
end