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
#############################################################################

function rameau(paths::Tuple,input_properties::Tuple,switches::Tuple;prediction_coef=0.069,temperature=23.0,laser=532.0,lb_break=2010.0,hb_start=1000.0,basetype="gcvspline")
	
	scale = 1000.
	
	liste=readcsv(paths[1], skipstart=1)
	names = liste[:,1]
	water = Array{Float64}(liste[:,3])
	smo_lf = Array{Float64}(liste[:,4])
	smo_hf = Array{Float64}(liste[:,5])
	roi = ones(Int((size(liste,2)-5)/2),2,size(liste,1))
	for j = 1:size(liste,1)
	    roi[:,1,j] = (liste[j,6:2:end])
	    roi[:,2,j] = (liste[j,7:2:end])
	end

	rws = ones(size(liste,1),3)
	
	for i = 1:size(liste,1)   
	    spectra = readdlm(string(paths[2],names[i]),input_properties[1],skipstart=input_properties[2])
    
	    if spectra[end,1] < spectra[1,1]
	        spectra = flipdim(spectra,1)
	    end
    
	    x = spectra[:,1]
	    y = spectra[:,2]./trapz(x,spectra[:,2]).*scale # area normalisation
		
		#### PRELIMINARY STEP: FIRST WE GRAB THE GOOD SIGNAL IN THE ROI
	    interest_index::Array{Int64} = find(roi[1,1,i] .<= x[:,1] .<= roi[1,2,i])
	    if size(roi)[1] > 1
	        for j = 2:size(roi)[1]
	            interest_index = vcat(interest_index,  find(roi[j,1,i] .<= x[:,1] .<= roi[j,2,i]))
	        end
	    end
	    interest_x = x[interest_index,1]
	    interest_y = y[interest_index,1]
    
	    if switches[2] == "yes" # if user asks for a double baseline
        
	        lb = roi[:,1,i]
	        hb = roi[:,2,i]

	        roi_hf =  [lb[hb.>hb_start] hb[hb.>hb_start]]
	        roi_lf =  [lb[hb.<lb_break] hb[hb.<lb_break]]
			
			#roi_lf = [lb[1,1] hb[1,1];x[y .== minimum(y[550 .< x .<900])]-5 x[y .== minimum(y[550 .< x .<900])]+5; lb[4,1] hb[4,1]]
        
	        # A first linear baseline in signals between ~1300 and ~ 2000 cm-1
	        y_calc1_part1, baseline1 = baseline(x,y,roi_hf[1,:],"poly",p=1.0)
			
			# we get a mean signal for background below 1300 cm-1
			X_low_lf = x[y .==minimum(y[1100.<x.<1400])]
			b_low_lf = mean(y[X_low_lf[1]-5 .< x .< X_low_lf[1]+5])
			baseline1[x.<X_low_lf[1]] = b_low_lf
			y_calc1_part1 = y - baseline1
        
	        # above 2500 cm-1, we use a thrid order polynomial function below the water peak
	        y_calc1_part2, baseline2 = baseline(x,y,roi_hf[2:end,:],basetype,p=[smo_hf[i]])
        	#y_calc1_part2, baseline2 = baseline(x,y,roi_hf[2:end,:],"poly",p=3.0)
		
	        # and we glue together the two parts
	        y_calc1 = [y_calc1_part1[0 .<x.<roi_hf[2,1]];y_calc1_part2[x .>= roi_hf[2,1]]]
	        bas1 = [baseline1[0 .<x.<roi_hf[2,1]];baseline2[x .>= roi_hf[2,1]]]
        
	        #y_calc1, bas1 = baseline(x,y,roi_hf,basetype=basetype,p=[smooth*10]) # baseline substraction
	        y_calc1 = y_calc1./trapz(x,y_calc1).*scale # area normalisation
        
	        if switches[3] == "yes"
	            x, y_long, ~ = tlcorrection([x[:] y_calc1[:]],temperature,laser) # Long correction
            	y_long = y_long .*scale
			
	            # Again we need a trick as we don't want to do anything in the HF part...
	            y_calc2_part1, bas2 = baseline(x,y_long,roi_lf,basetype,p=[smo_lf[i]]) # baseline substraction
				#y_calc2_part1, bas2 = baseline(x,y_long,roi_lf,"poly",p=4.0) # baseline substraction
	            y_calc2 = [y_calc2_part1[x.<roi_lf[end,1]];y_long[x.>=roi_lf[end,1]]] # we glue the good parts 
	            bas2[x.>roi_lf[end,1]] = 0.0 # we put this part to 0
	            y_calc2 = y_calc2./trapz(x,y_calc2).*scale # area normalisation
	        else
            
	            y_calc2_part1, bas2 = baseline(x,y_calc1,roi_lf,basetype,p=[smo_lf[i]]) # baseline substraction
	            y_calc2 = [y_calc2_part1[x.<roi_lf[end,1]];y_calc1[x.>=roi_lf[end,1]]] # we glue the good parts
	            bas2[x.>roi_lf[end,1]] = 0.0 # we put this part to 0
	            y_calc2 = y_calc2./trapz(x,y_calc2).*scale # area normalisation
            
	        end
			
	        # Plotting the spectra
	        figure(figsize=(10,10))
			suptitle("Spectrum $(names[i])")
			
	        subplot(221)
	        plot(x,y,color="black",label="Raw sp.")
			scatter(interest_x,interest_y,s=10.,color="red",label="ROI")
	        plot(x,bas1,color="blue",label="Baseline")
	        xlabel(L"Raman shift, cm$^{-1}$")
	        ylabel("Intensity, area normalised")
	        title("First baseline")
			legend(loc="best",fancybox="true")
			
	        subplot(223)
	        plot(x,y,color="black",label="Raw sp.")
			scatter(interest_x,interest_y,s=10.,color="red",label="ROI")
	        plot(x,bas1,color="blue",label="Baseline")
			xlim(2300,4000)
			ylim(minimum(y[x.>2300]),maximum(y[x.>2300]))
	        xlabel(L"Raman shift, cm$^{-1}$")
	        ylabel("Intensity, area normalised")
	        title("First baseline, HF zoom")
        
	        subplot(222)
	        if switches[3] == "yes" # if Long correction is activated
	            plot(x,y_long,color="blue",label="Long corr. sp.")
				scatter(interest_x,y_long[interest_index,1],s=10.,color="red",label="ROI") # the ROI signals after first baseline and long correction
	        else
	            plot(x,y_calc1,color="blue",label="Corr. sp.")
				scatter(interest_x,y_calc1[interest_index,1],s=10.,color="red",label="ROI") # the ROI signals after first baseline
	        end
	        plot(x,bas2,color="magenta",label="2nd baseline")
	        plot(x,y_calc2,color="cyan",label="Final sp.")
	        xlabel(L"Raman shift, cm$^{-1}$")
	        ylabel("Intensity, area normalised")
	        title("Second baseline")
			legend(loc="best",fancybox="true")
			
	        subplot(224)
	        if switches[3] == "yes" # if Long correction is activated
	            plot(x,y_long,color="blue",label="Long corr. sp.")
				scatter(interest_x,y_long[interest_index,1],s=10.,color="red",label="ROI") # the ROI signals after first baseline and long correction
				ylim(0,maximum(y_long[x.<1500]))
	        else
	            plot(x,y_calc1,color="blue",label="Corr. sp.")
				scatter(interest_x,y_calc1[interest_index,1],s=10.,color="red",label="ROI") # the ROI signals after first baseline
				ylim(0,maximum(y_calc1[x.<1500]))
	        end
	        plot(x,bas2,color="magenta",label="2nd baseline")
	        plot(x,y_calc2,color="cyan",label="Final sp.")
			xlim(0,1500)
			
	        xlabel(L"Raman shift, cm$^{-1}$")
	        ylabel("Intensity, area normalised")
	        title("Second baseline, LF zoom")
			
            
	    else # only a single baseline treatment is asked
        
	        # we test if thre Long correction is asked, and react accordingly for baseline substraction
	        if switches[3] == "yes"
	            x, y_long, ~ = tlcorrection([x[:] y[:]],temperature,laser) # Long correction
            	#y_long = y_long ./maximum(y_long)
				y_long=y_long.*scale
				
	            y_calc2, bas2 = baseline(x,y_long,roi[:,:,i],basetype,p=[smo_lf[i]])
				y_calc2 = y_calc2[:,1]./trapz(x[:],y_calc2[:,1]).*(scale./10) # area normalisation
            
	            figure()
	            plot(x,y,"black",label="Raw sp.")
	            plot(x,y_long,"blue",label="Long corr. sp.")
				scatter(interest_x,y_long[interest_index,1],s=10.,color="red",label="ROI") # the ROI signals after first baseline and long correction
	            plot(x,bas2,"red",label="Baseline")
	            plot(x,y_calc2,"cyan",label="Final sp.")
	            xlabel(L"Raman shift, cm$^{-1}$")
	            ylabel("Intensity, area normalised")
	            title("Baseline sub. with Long corr., sp. $(names[i])")
				legend(loc="best",fancybox="true")
            
	        else
	            y_calc2, bas2 = baseline(x,y,roi[:,:,i],basetype,p=[smo_lf[i]])
				y_calc2 = y_calc2[:,1]./trapz(x[:],y_calc2[:,1]).*scale # area normalisation
            
	            figure()
	            plot(x,y.*10,"black",label="Raw sp.")
				scatter(interest_x,interest_y.*10,s=10.,color="red",label="ROI")
	            plot(x,bas2.*10,"red",label="Baseline")
	            plot(x,y_calc2,"cyan",label="Final sp.")
	            xlabel(L"Raman shift, cm$^{-1}$")
	            ylabel("Intensity, area normalised")
	            title("Baseline sub., NO Long corr.")
				legend(loc="best",fancybox="true")
	        end
        
	    end
    
	    # saving the corrected spectrum
	    writecsv(string(paths[3],names[i]),[x y_calc2 bas2])
    
	    #Saving the generated figure (specified directory should exist)
	    savefig(string(paths[4],"baseline_",names[i],".pdf"))
    
	    # calc=ulating the areas under silicate and water bands
	    As = trapz(x[200 .< x .<1300],y_calc2[200 .< x .<1300])
	    Aw = trapz(x[3100 .< x .<3750],y_calc2[3100 .< x .<3750])
    
	    # recording them
	    rws[i,1] = As;
	    rws[i,2] = Aw;
	    rws[i,3] = Aw./(As);
	    println("The rws is $(rws[i,3])")
	end
	
	if switches[1] == "yes"
	    calibration_model(x, p) = p[1].*x
	    fit = curve_fit(calibration_model, rws[:,3], water[:]./(100.-water[:]), [0.5])
	    coef = fit.param

	    rws_calibration = collect(0:0.01:round(maximum(rws[:,3]),2))
	    water_ratio_calibration = calibration_model(rws_calibration,coef)
	    water_compare =  100.*calibration_model(rws[:,3],coef)./(calibration_model(rws[:,3],coef)+1) # eq. 3 Le Losq et al. (2012)

	    rmse_calibration = sqrt(1./(size(rws,1)-1).*sum((water_compare-water).^2))
        
	    figure()
	    scatter(rws[:,3],water[:]./(100.-water[:]))
	    plot(rws_calibration,water_ratio_calibration,color="red")
	    xlabel("Awater/Asilicates, area ratio")
	    ylabel("Water/Glass, weight ratio")
	    title("Calibration with a Rws factor of $(round(coef,5)), water% std = $(round(rmse_calibration,2))")
	    savefig(paths[6])
        
	    writecsv(string(paths[5]), [rws[:,3] water water_compare]   )
	    #return [rws[:,3] water water_compare]   
	else
	    water_predicted = ws[:,3].*provided_coef 
	    writecsv(string(paths[5]), [rws[:,3] water_predicted])
	    #return [rws[:,3] water_predicted]
	end           
end

