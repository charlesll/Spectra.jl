#############################################################################
#Copyright (c) 2016 Charles Le Losq
#
#The MIT License (MIT)
#
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the #Software without restriction, including without limitation the rights to use, copy, #modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, #and to permit persons to whom the Software is furnished to do so, subject to the #following conditions:
#
#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, #INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
# diffusion.jl contains several functions focused on the treatment of infrared spectra along diffusion profiles in minerals
#
#
#############################################################################
"""
The 'IRdataprep' function allows to obtain the frequency, distance and y vectors in a region of interest, ready to input in an optimisation algorithm.

   IRdataprep(data,distance_step,start_sp,stop_sp,low_x,high_x,norm_low_x,norm_high_x)
with
data (array, Float64): an array containing in the first column the x values and in the subsequent columns the correcponding intensities;
distance_step (Float64): the steps between each spectrum along the diffusion profile, in metres;
start_sp and stop_sp (Int): the starting and stopping spectra of interest (warning: number of data columns - 1, as data also contains the x axis such that size(data)[2] = number of spectras + 1);
low_x,high_x (Float64): frequencies between which the signal is of interest, defined as low_x < high_x;
norm_low_x and norm_high_x (Float64): frequencies between which the signal must be integrated for mormalising the final signals.

This function first checks that the x values are increasing, and correct that if this is not the case.

Then, it selects the data in the region of interest, delimited by frequencies defined as low_x < high_x.
Normalisation to the area between the norm_low_x and norm_high_x frequency is also performed (norm_low_x < norm_high_x)

The code returns the x and the y arrays, plus arrays for optimisation: x_input_fit contains the frequency-distance couples and y_input_fit their corresponding y values.
"""

function IRdataprep(data::Array{Float64},distance_step::Float64,start_sp::Int,stop_sp::Int,low_x::Float64,high_x::Float64,norm_low_x::Float64,norm_high_x::Float64)
   start_sp = start_sp+1
   stop_sp = stop_sp+1
   
   if data[end,1] <= data[1,1]
      data=data[end:-1:1,:]
   end
   
   x::Array{Float64} = data[find(low_x .< data[:,1] .< high_x),1]
   x_sili::Array{Float64} = data[find(norm_low_x .< data[:,1] .< norm_high_x),1]

   if (start_sp+stop_sp) .> size(x)[1]
      y::Array{Float64} = data[find(low_x .< data[:,1] .< high_x),start_sp:end]
      y_sili::Array{Float64} = data[find(norm_low_x .< data[:,1] .< norm_high_x),start_sp:end]
   else
      y = data[find(low_x .< data[:,1] .< high_x),start_sp:start_sp+stop_sp]
      y_sili = data[find(norm_low_x .< data[:,1] .< norm_high_x),start_sp:start_sp+stop_sp]
   end

   # constructing the good distance vector + x and Y associated values
   steps::Float64 = 0.0
   x_input_fit::Array{Float64} = zeros(length(x),2)
   y_input_fit::Array{Float64} = zeros(length(x),1)
   for i = 1:size(y)[2]
      x_temp::Array{Float64} = zeros(length(x),2)
      if i == 1
         x_input_fit[:,1] = steps
         x_input_fit[:,2] = x
         y[:,i] = y[:,i]/trapz(x_sili[:,1],y_sili[:,i]) # integration of silicate bands for normalisatio
         y_input_fit[:,1] = y[:,i]
      else
         x_temp[:,1] = steps
         x_temp[:,2] = x
         y[:,i] = y[:,i]/trapz(x_sili[:,1],y_sili[:,i]) # integration of silicate bands for normalisatio
        
         x_input_fit = vcat(x_input_fit,x_temp)
         y_input_fit = vcat(y_input_fit,y[:,i])
   end
   steps = steps + distance_step
   end
   
   return x, y, x_input_fit, y_input_fit
end

"""
The `peak_diffusion` function allows to define a gaussian peak for which the amplitude is defined by a 1D diffusion law
Call as:
    peak_diffusion(distance,temps,w,c0,c1,diffusion,w0,hwhm)

where:
C1 is the concentration imposed, at the diffusive bondary (a.u., float64 number);
C0 is the initial concentration in the host (a.u., float64 number);
Distance is the distance in metre from the diffusive bondary (Vector float 64);
temps is the diffusive time in seconds (float64 number);
w is the measured frequency (Vector float 64);
w0 is the peak centre frequency (float64 number);
hwhm is the peak half-width at middle height (float64 number).

The function returns a vector of the values of the peak intensity at the input frequencies.

The current version assumes the peak being a Gaussian peak. Updates will integrate Lorentzian and Pseudo-Voigt peaks.
"""
function peak_diffusion(distance::Array{Float64,1},temps::Float64,w::Array{Float64,1},c0::Float64,c1::Float64,diffusion::Float64,w0::Float64,hwhm::Float64)
    amplitude::Array{Float64,1} = (c1 - c0) * erfc(distance[:,1] ./ (2. * sqrt( 10^diffusion * temps))) + c0
    gaussienne::Array{Float64,1} = exp(-log(2) .* ( (w[:,1]-w0) ./ hwhm ).^2)
    diff_gauss =amplitude .* gaussienne
    return diff_gauss
end

"""
The global model for using with the Julia.optim package. Not necessarily up to date because I prefer using the JuMP way.
   model(p,x_fit,temps;compo_out = 0)
"""
function model(p::Array{Float64,1},x_fit::Array{Float64,2},temps::Float64;compo_out::Int = 0)
    number_gauss = length(p)/5
    components = zeros(length(x_fit[:,2]),round(Int8,number_gauss))
    compteur = 1
    for i = 1:5:length(p)
        components[:,compteur] = peak_diffusion(x_fit[:,1],temps,x_fit[:,2],p[i],p[i+1],p[i+2],p[i+3],p[i+4])
        compteur = compteur + 1
    end
    if compo_out == 0
        return sum(components,2)
    else
        return components
    end
end

