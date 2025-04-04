#############################################################################
#Copyright (c) 2016-2025 Charles Le Losq
#
#The MIT License (MIT)
#
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the #Software without restriction, including without limitation the rights to use, copy, #modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, #and to permit persons to whom the Software is furnished to do so, subject to the #following conditions:
#
#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, #INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
# functions.jl contains several mathematic functions
#
#
#############################################################################

"""
	bootsample(x::Array{Float64}, y::Array{Float64}; boottype::String = "np", ese::Array{Float64} = [0.0])

Inputs
------

	x: Array{Float64}, the x axis. It can have multiple columns.

	y: Array{Float64}, the y axis. It can have multiple columns.

Options
-------

	boottype: ASCIIString, either "np" or "p", this is the type of bootstrapping performed. "np" performes a non-parametric resampling fo the dataset with replacement. "p" performs a parametric resampling. The data are resample from a gaussian distribution centered on the y values with errors that should be provided in the ese variable.

	ese: Array{Float64}, containing the errors affecting the y values that are used during parametric bootstrapping.

Outputs
-------

	b_x_f: Array{Float64}, the bootstrapped x values

	b_y_f: Array{Float64}, the bootstrapped y values

The bootstrap function can be embedded in a for loop, and will each time produce a different dataset. Performing K times the bootstrapping and fitting each time the model will allow to estimate the error distribution on the peak parameters. This technic has the advantage of making no prior assumption on the probability distribution functions of parameters errors. However, it is  much more time consuming that using the covariance matrix.

"""
function bootsample(x::Array{Float64}, y::Array{Float64}; boottype::String = "np", ese::Array{Float64} = [0.0])
    testx = size(x,1)
    testy = size(y,1)
    teste = size(ese,1)
    if boottype == "np" && testx == testy
        vect = collect(1:size(x)[1]) # for real bootstrapping
        idx = sample(vect,size(vect)[1],replace=true) #resampling data with replacement...
        b_x_f = x[idx,:] # we pick the right x
        b_y_f = y[idx,:] # we pick the right y
    elseif boottype == "p" && testx == testy && testx == teste
	    b_x_f = x
	    b_y_f = ones(size(y))
	    b_y_f = randn!(b_y_f[:,:]) .* ese[:,:] + y[:,:]
    else
	    error("Something is wrong. Check size of x, y and ese as well as the boottype (either p or np). Providing an ese array is mandatory for parametric bootstrapping")
    end

    return b_x_f, b_y_f
end

"""
	bootperf(params_boot::Array{Float64}; plotting::String = "True", parameter::Int64 = 1, feature::Int64 = 1, histogram_step::Int64 = 100, savefigures::String = "False", save_bootrecord::String = "Boot_record.pdf", save_histogram::String = "Boot_histogram.pdf")

Inputs
------

	params_boot[i,j] or [i,j,k]
		array with i the bootstrrap experiment, j the parameter and k the feature being calculated (e.g., a peak during peak fitting)
	plotting: String
		switch to plotting mode ("True") or not ("False"). If true, parameter and feature must be provided, otherwise an error message is returned.

	histogram_step
		integer value to control the histogram X axis division.
	savefigures: "True" or "False"
		explicit, save in the current working directory.
	save_bootrecord:
		Name for the graphic showing the bootstrap mean and std evolutions, with extension.
	save_histogram:
		Name for the graphic showing the histogram for the parameter and feature of interest, with extension.

Outputs
-------

	std_record, mean_record, the arrays recording how the standard deviation and mean of the parameters as a function of the bootstrap advance.

"""
function bootperf(params_boot::Array{Float64}; plotting::String = "True", parameter::Int64 = 1, feature::Int64 = 1, histogram_step::Int64 = 100, savefigures::String = "False", save_bootrecord::String = "Boot_record.pdf", save_histogram::String = "Boot_histogram.pdf")

    # test
    if plotting != "True" && plotting != "False"
        error("Error: plotting should be set to True or False")
    end

    nb_boot::Int64 = size(params_boot)[1]
    nb_param::Int64 = size(params_boot)[2]

   if ndims(params_boot) == 2 #### TWO DIMENSIONS PARAMS_BOOT

       mean_record::Array{Float64} = zeros(size(params_boot))
       std_record::Array{Float64} = zeros(size(params_boot))
       for k = 1:nb_boot
           for j = 1:size(params_boot)[2]
       	        if k == 1
       	            mean_record[k,j] = params_boot[k,j]
       		    std_record[k,j] = 0
       	        else
       	            mean_record[k,j] = sum(mean(params_boot[1:k,j],1))
       		    std_record[k,j] = sum(std(params_boot[1:k,j],1))
       	        end
       	    end
       	end


       if plotting == "True" # Graphics
           # tests:
           if savefigures != "True" && savefigures != "False"
               error("savefigure should be set to True or False.")
           end
           if isa(save_bootrecord ,String) != true
               error("Error: save_bootrecord is not a valid savename. Check it is a string of characters.")
           end
           if isa(save_histogram ,String) != true
               error("Error: save_histogram is not a valid savename. Check it is a string of characters.")
           end

           if parameter == 0 || histogram_step == 0
               error("parameter or histogram_step appear to have been set to 0...")
           elseif parameter > size(params_boot)[2]
               error("Value entered for parameter is too high.")
           end


           fig = figure()
           #title("Mean and Std of parameter $(parameter), feature $(feature)",fontsize = 10, fontweight = "bold", fontname="Arial")
           subplot(121)
           plot(1:nb_boot, mean_record[:,parameter],label="Mean")
           xlabel("Number of iterations during bootstrap",fontsize = 10, fontweight = "bold", fontname="Arial")
           ylabel("Mean value",fontsize = 10, fontweight = "bold", fontname="Arial")

           subplot(122)
           plot(1:nb_boot, std_record[:,parameter],label="Std")
           xlabel("Number of iterations during bootstrap",fontsize = 10, fontweight = "bold", fontname="Arial")
           ylabel("Standard deviation",fontsize = 10, fontweight = "bold", fontname="Arial")

           if savefigures == "True"
               savefig(save_bootrecord)
           end

           fig2= figure()
           h = plt[:hist](params_boot[:,parameter],100)
           #title("Histogram of parameter $(parameter), feature $(feature)",fontsize = 10, fontweight = "bold", fontname="Arial")
           xlabel("Parameter value",fontsize = 10, fontweight = "bold", fontname="Arial")
           ylabel("Number of bootstrap values",fontsize = 10, fontweight = "bold", fontname="Arial")

           if savefigures == "True"
               savefig(save_histogram)
       	   end

       end # end if plotting

       return mean_record, std_record


   elseif ndims(params_boot) == 3 #### THREE DIMENSIONS PARAMS_BOOT
       mean_record = zeros(size(params_boot))
       std_record = zeros(size(params_boot))
       for k = 1:nb_boot
           for j = 1:size(params_boot)[2]
   	    for i = 1:size(params_boot)[3]
   	        if k == 1
   	            mean_record[k,j,i] = params_boot[k,j,i]
   		    std_record[k,j,i] = 0
   	        else
   	            mean_record[k,j,i] = sum(mean(params_boot[1:k,j,i],1))
   		    std_record[k,j,i] = sum(std(params_boot[1:k,j,i],1))
   	        end
   	    end
   	end
       end

       if plotting == "True" # Graphics
           # tests:
           if savefigures != "True" && savefigures != "False"
               error("savefigure should be set to True or False.")
           end
           if isa(save_bootrecord ,String) != true
               error("Error: save_bootrecord is not a valid savename. Check it is a string of characters.")
           end
           if isa(save_bootrecord ,String) != true
               error("Error: save_histogram is not a valid savename. Check it is a string of characters.")
           end

           if parameter == 0 || feature == 0 || histogram_step == 0
               error("parameter, feature or histogram_step appear to have been set to 0...")
           elseif parameter > size(params_boot)[2]
               error("Value entered for parameter is too high.")
           elseif feature > size(params_boot)[3]
   	    error("Value entered for feature is too high.")
           end


           fig = figure()
           #title("Mean and Std of parameter $(parameter), feature $(feature)",fontsize = 10, fontweight = "bold", fontname="Arial")
           subplot(121)
           plot(1:nb_boot, mean_record[:,parameter,feature],label="Mean")
           xlabel("Number of iterations during bootstrap",fontsize = 10, fontweight = "bold", fontname="Arial")
           ylabel("Mean value",fontsize = 10, fontweight = "bold", fontname="Arial")

           subplot(122)
           plot(1:nb_boot, std_record[:,parameter,feature],label="Std")
           xlabel("Number of iterations during bootstrap",fontsize = 10, fontweight = "bold", fontname="Arial")
           ylabel("Standard deviation",fontsize = 10, fontweight = "bold", fontname="Arial")

           if savefigures == "True"
               savefig(save_bootrecord)
           end

           fig2= figure()
           h = plt[:hist](params_boot[:,parameter,feature],100)
           #title("Histogram of parameter $(parameter), feature $(feature)",fontsize = 10, fontweight = "bold", fontname="Arial")
           xlabel("Parameter value",fontsize = 10, fontweight = "bold", fontname="Arial")
           ylabel("Number of bootstrap values",fontsize = 10, fontweight = "bold", fontname="Arial")

           if savefigures == "True"
               savefig(save_histogram)
   	end

       end # end if plotting

       return mean_record, std_record

    elseif ndims(params_boot) > 3
	error("Not implemented, params_boot should have 3 dimensions, maximum.")
    end

end # end function
