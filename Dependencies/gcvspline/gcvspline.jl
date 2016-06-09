using PyCall
unshift!(PyVector(pyimport("sys")["path"]), "")

@pyimport gcvspline

function Gspline(xfit::Array{Float64},yfit::Array{Float64},esefit::Array{Float64},xtarget::Array{Float64},smoothing::Array{Float64};VAL::Array{Float64}= esefit.^2,splorder::Float64=2.,splmode::Float64=3.,NC::Float64 = size(yfit,1))
	
    c, wk, ier = gcvspline.gcvspline(xfit,yfit,smoothing.*ese,VAL,splmode, splorder,NC) # gcvspl with mode 3 and smooth factor
    out::Array{Float64} = gcvspline.splderivative(x,xfit,c,splorder; L = 1,IDER = 0) 
	
end