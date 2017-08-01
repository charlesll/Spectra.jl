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

"""
	rameau(paths::Tuple,switches::Tuple;input_properties=('\t',0),prediction_coef=[0.0059;0.0005],temperature=23.0,laser=532.0,lb_break=1600.,hb_start=2600.,roi_hf_external = [3000. 3100.; 3800. 3900.],basetype="gcvspline",mmap_switch=true)
	
INPUTS:

    paths: Tuple{Strings}, it contains the following strings:

		in_liste: the relative path of your liste of spectra (see the description of this liste below), e.g. "./liste_2012.csv"

		in_path = the relative path of the repertory where your raw spectra are stored, e.g. "./raw/"

		out_path = the relative path of the repertory where you want to output the corrected spectra, e.g. "./treated/" (THIS FOLDER MUST EXIST PRIOR TO RUN)

		fig_path= the relative path of the repertory for saving the figures, e.g. "./figures/"

		rws_save_file = the relative path of the csv file where the results will be stored, e.g. "rws_calculated.csv". If calibration is set to "yes", it outputs an array containing the Rws, the provided water content and the calculated water content in columns 1, 2 and 3, respectively. If calibration is set to anything else, it outputs the rws and the water content predicted with the provided prediction_coef

		rws_save_fig = the relative path of the figure for the calibration, e.g. "rws_calculated.jpg". Note: It only works when calibration is on "yes", and Type_of_calibration is set on "internal". Just set it to "trash.jpg" if you are not going to use it.

	input_properties: Tuple, it contains the delimiter of the columns in your spectra files, and the number of starting lines to skip. For instance, use ('\t',1) if your files have one header line and the columns are separated  by tabulations, or (',',0) if you use CSV files without header.

	switches: Tuple, it contains the different switches as (Type_of_calibration, calibration?, experimental?, temperature_laser_correction?).

		Type_of_calibration: should be set to "internal" or "external". The external mode requires standards, as described in Behrens et al. (2006) and Thomas et al. (2008). If you use the external mode, you can leave blanks for the other switches. For instance, enter ("external","","","").

		calibration? :  use only with the "internal" mode, this should be equal to "yes" or "no". Enter yes if you are setting up a new calibration, or "no" if you want to use a previously determined calibration coefficient with providing the later in the prediction_coef variable (see options).

		experimental? : use only with the "internal" mode, this is an experimental code. For now, only the "double" feature is advised. It allows using different smoothing spline coefficients for the water and silicate bands, delimited by the lb_break and hb_start variables (see below).

		temperature_laser_correction? : use only with the "internal" mode, this should be equal to "yes" or "no". This asks if you want to use the temperature-laser wavelength correction as done in Le Losq et al. (2012). If you use the double baseline mode, this correction is applied after removing the background under the water peak.

OPTIONS:

	prediction_coef: Array{Float64}, this array contains the calibration coefficient with its error bar. Those values will be used in the predictions, if you use the predictive mode (i.g., calibration switch is not set to "yes"). Default = [0.0059;0.0005].

	temperature: Float64, the temperature for the temperature-laser wavelength correction in Celsius. Default = 23.0.

	laser: Float64, the laser wavelength in nm for the temperature-laser wavelength correction. Default = 532.0.

	lb_break: Float64, for double baseline correction, the breaking point before which the software will consider the BIRs in the low frequency region. Default = 2010.0.

	hb_start: Float64, for double baseline correction, the breaking point after which the software will consider the BIRs in the high frequency region. Default = 1000.0.

  roi_hf_external: Array{Float64}, the roi for fitting the linear baseline in the external calibration mode. Default = [3000. 3100.; 3800. 3900.].

	basetype: String, the type of baseline you want to fit. Corresponds to the "basetype" parameter of the baseline function.  Default = "gcvspline".

    mmap_switch: false or true, this allows to switch on or off the memory mapping in the `readcsv`/`readdlm` functions that `rameau` uses. Default = "true".

OUTPUTS:

	Rameau does not provide any outputs directly in Julia, but saves everything in the folders you indicate in the variable "paths".

"""


function rameau(paths::Tuple,switches::Tuple;input_properties=('\t',0),prediction_coef=[0.0059;0.0005],temperature=23.0,laser=532.0,lb_break=1600.,hb_start=2600.,roi_hf_external = [3000. 3100.; 3800. 3900.],basetype="gcvspline",mmap_switch=true,lambda=10.0^12,p=0.01)

	# some function definition
	calibration_model(x, p) = p[1].*x
	baseline_lfmodel(x, p) = p[1].*x + p[2].*x.^2 + p[3].*x.^3

	if switches[1] == "internal"


		scale = 1 # no scaling, keeping that for history

		liste=readcsv(paths[1], skipstart=1,use_mmap = mmap_switch)
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
			spectra = Array{Float64}(readdlm(string(paths[2],names[i]),input_properties[1],skipstart=input_properties[2],use_mmap = mmap_switch))

			if spectra[end,1] < spectra[1,1]
				spectra = flipdim(spectra,1)
			end

			x = spectra[:,1]
			spectra[:,2] = spectra[:,2] - minimum(spectra[:,2])
			y = spectra[:,2]./maximum(spectra[:,2]).*scale + 0.1# area normalisation, and we avoid any 0 by putting the smallest value to 1e-20
			ratio_bkg = minimum(y)/maximum(y) # to keep a record of the ratio of maximum signal intensity over minimum background intensity

			#### PRELIMINARY STEP: FIRST WE GRAB THE GOOD SIGNAL IN THE ROI
			interest_index::Array{Int64} = find(roi[1,1,i] .<= x[:,1] .<= roi[1,2,i])
			if size(roi)[1] > 1
				for j = 2:size(roi)[1]
					interest_index = vcat(interest_index,  find(roi[j,1,i] .<= x[:,1] .<= roi[j,2,i]))
				end
			end
			interest_x = x[interest_index,1]
			interest_y = y[interest_index,1]

			if switches[3] == "double" # if user asks for the experimental mode

				if switches[4] == "yes"

					# If the Long correction is asked, then it will use the two smoothing spline factors and perform a "traditional" internal baseline fit with gcvspline
					x, y_long, ~ = tlcorrection([x[:] y[:]],temperature,laser) # Long correction
					y_long=y_long./maximum(y_long).*scale

					lb = roi[:,1,i]
					hb = roi[:,2,i]

					roi_hf =  [lb[hb.>hb_start] hb[hb.>hb_start]]
					roi_lf =  [lb[hb.<lb_break] hb[hb.<lb_break]]

					y_hf = y_long[x.>hb_start]
					y_lf = y_long[x.<lb_break]

					# A first linear baseline in signals between ~1300 and ~ 2000 cm-1
					y_calc_hf, baseline_hf = baseline(x,y_long,roi_hf,basetype,p=smo_hf[i])
					y_calc_lf, baseline_lf = baseline(x,y_long,roi_lf,basetype,p=smo_lf[i])

					y_calc2 = [y_calc_lf[x.<=lb_break];y_calc_hf[x.>lb_break]]
					bas2 = [baseline_lf[x.<=lb_break];baseline_hf[x.>lb_break]]
					#y_calc2 = y_calc2[:,1]./trapz(x[:],y_calc2[:,1]).*scale # area normalisation

					figure(figsize=(20,20))
					suptitle("LL2012 method, sp. $(names[i])")

					subplot(3,2,(1,2))
					plot(x,y,"black",label="Raw sp.")
					plot(x,y_long,"blue",label="Long corr. sp.")
					scatter(interest_x,y_long[interest_index,1],s=10.,color="red",label="ROI") # the ROI signals after first baseline and long correction
					plot(x,bas2,"red",label="Baseline")
					legend(loc="best",fancybox="true")

					xlabel(L"Raman shift, cm$^{-1}$",fontsize=18,fontname="Arial")
					ylabel("Intensity,a. u.",fontsize=18,fontname="Arial")

					subplot(3,2,3)
					plot(x[x.<1500],y_long[x.<1500],"blue",label="Long corr. sp.")
					scatter(interest_x[interest_x.<1500],y_long[interest_index[interest_x.<1500],1],s=10.,color="red",label="ROI") # the ROI signals after first baseline and long correction
					plot(x[x.<1500],bas2[x.<1500],"red",label="Baseline")
					xlim(0,1500)
					legend(loc="best",fancybox="true")
					xlabel(L"Raman shift, cm$^{-1}$",fontsize=18,fontname="Arial")
					ylabel("Intensity,a. u.",fontsize=18,fontname="Arial")

					subplot(3,2,4)
					plot(x[x.>2800],y_long[x.>2800],"blue",label="Long corr. sp.")
					scatter(interest_x[interest_x.>2800],y_long[interest_index[interest_x.>2800],1],s=10.,color="red",label="ROI") # the ROI signals after first baseline and long correction
					plot(x[x.>2800],bas2[x.>2800],"red",label="Baseline")
					xlim(2600,4000)
					legend(loc="best",fancybox="true")
					xlabel(L"Raman shift, cm$^{-1}$",fontsize=18,fontname="Arial")
					ylabel("Intensity,a. u.",fontsize=18,fontname="Arial")

					subplot(3,2,(5,6))
					plot(x,y_calc2,"cyan",label="Final sp.")
					xlabel(L"Raman shift, cm$^{-1}$")
					ylabel("Intensity, area normalised")
					legend(loc="best",fancybox="true")
					xlabel(L"Raman shift, cm$^{-1}$",fontsize=18,fontname="Arial")
					ylabel("Intensity,a. u.",fontsize=18,fontname="Arial")

				else
					# ALSO FULLY EXPERIMENTAL!!!! IT DOES NOT WORK!
					y_hf = y[x.>hb_start]
					y_lf = y[x.<lb_break]

					# A first linear baseline in signals between ~1300 and ~ 2000 cm-1
					y_calc_hf, baseline_hf = baseline(x,y,[2600. 2900.;3775. 3900.],"poly",p=3.0)

				    #A first linear baseline in signals between ~1300 and ~ 2000 cm-1
				   	y_calc1_part1, baseline1 = baseline(x,y,[1300 1800.],"poly",p=1.0)

					# we get a mean signal for background below 1300 cm-1
					X_low_lf = x[y .== minimum(y[100 .<x.<1400])]
					b_low_lf = mean(y[X_low_lf[1]-10 .< x .< X_low_lf[1]+10])

					baseline1[x.<1250.] = b_low_lf
					y_calc1_part1 = y - baseline1

					bas2 = [baseline1[x.<lb_break];baseline_hf[x.>lb_break]]
					y_calc2 = [y_calc1_part1[x.<lb_break];y_calc_hf[x.>lb_break]]

					figure()
					plot(x,y,color="black")
					plot(x,baseline_hf,color="red")
					plot(x,y_calc2,color="blue")
				end

			else # only a single baseline treatment is asked

				# we test if the Long correction is asked, and react accordingly for baseline substraction
				if switches[4] == "yes"
					x, y_long, ~ = tlcorrection([x[:] y[:]],temperature,laser) # Long correction
					y_long=y_long./maximum(y_long).*scale
					
					if basetype == "arPLS" || basetype == "whittaker"
						y_calc2, bas2 = baseline(x,y_long,roi[:,:,i],basetype,lambda = lambda, p = p)
					else
						y_calc2, bas2 = baseline(x,y_long,roi[:,:,i],basetype,p=smo_lf[i])
					end
					
					y_calc2 = y_calc2[:,1]./trapz(x[:],y_calc2[:,1]).*(scale./10) # area normalisation

					figure(figsize=(20,20))
					suptitle("LL2012 method, sp. $(names[i])")

					subplot(3,2,(1,2))
					plot(x,y,"black",label="Raw sp.")
					plot(x,y_long,"blue",label="Long corr. sp.")
					scatter(interest_x,y_long[interest_index,1],s=10.,color="red",label="ROI") # the ROI signals after first baseline and long correction
					plot(x,bas2,"red",label="Baseline")
					legend(loc="best",fancybox="true")

					xlabel(L"Raman shift, cm$^{-1}$",fontsize=18,fontname="Arial")
					ylabel("Intensity,a. u.",fontsize=18,fontname="Arial")

					subplot(3,2,3)
					plot(x[x.<1500],y_long[x.<1500],"blue",label="Long corr. sp.")
					scatter(interest_x[interest_x.<1500],y_long[interest_index[interest_x.<1500],1],s=10.,color="red",label="ROI") # the ROI signals after first baseline and long correction
					plot(x[x.<1500],bas2[x.<1500],"red",label="Baseline")
					xlim(0,1500)
					legend(loc="best",fancybox="true")
					xlabel(L"Raman shift, cm$^{-1}$",fontsize=18,fontname="Arial")
					ylabel("Intensity,a. u.",fontsize=18,fontname="Arial")

					subplot(3,2,4)
					plot(x[x.>2800],y_long[x.>2800],"blue",label="Long corr. sp.")
					scatter(interest_x[interest_x.>2800],y_long[interest_index[interest_x.>2800],1],s=10.,color="red",label="ROI") # the ROI signals after first baseline and long correction
					plot(x[x.>2800],bas2[x.>2800],"red",label="Baseline")
					xlim(2600,4000)
					legend(loc="best",fancybox="true")
					xlabel(L"Raman shift, cm$^{-1}$",fontsize=18,fontname="Arial")
					ylabel("Intensity,a. u.",fontsize=18,fontname="Arial")

					subplot(3,2,(5,6))
					plot(x,y_calc2,"cyan",label="Final sp.")
					xlabel(L"Raman shift, cm$^{-1}$")
					ylabel("Intensity, area normalised")
					legend(loc="best",fancybox="true")
					xlabel(L"Raman shift, cm$^{-1}$",fontsize=18,fontname="Arial")
					ylabel("Intensity,a. u.",fontsize=18,fontname="Arial")

				else
					if basetype == "arPLS" || basetype == "whittaker"
						y_calc2, bas2 = baseline(x,y_long,roi[:,:,i],basetype,lambda = lambda, p = p)
					else
						y_calc2, bas2 = baseline(x,y_long,roi[:,:,i],basetype,p=smo_lf[i])
					end
					
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
    		close()
			if switches[3] == "double"
				# calc=ulating the areas under silicate and water bands
				As = trapz(x[150 .< x .<1300],y_calc2[150 .< x .<1300])
				Aw = trapz(x[2800 .< x .<3800],y_calc2[2800 .< x .<3800])
			else
				# calc=ulating the areas under silicate and water bands
				As = trapz(x[150 .< x .<1300],y_calc2[150 .< x .<1300])
				Aw = trapz(x[2800 .< x .<3800],y_calc2[2800 .< x .<3800])
			end

			# recording them
			rws[i,1] = As;
			rws[i,2] = Aw;
			rws[i,3] = Aw./(As);
			println("Spectrum $(names[i]), the rws is $(rws[i,3])")
		end

		if switches[2] == "yes"
			fit = curve_fit(calibration_model, rws[:,3], water[:]./(100.-water[:]), [0.5])
			coef = fit.param
			sigma = estimate_errors(fit, 0.95)

			rws_calibration = collect(0:0.01:round(maximum(rws[:,3]),2))
			water_ratio_calibration = calibration_model(rws_calibration,coef)
			water_compare =  100.*calibration_model(rws[:,3],coef)./(calibration_model(rws[:,3],coef)+1) # eq. 3 Le Losq et al. (2012)

			rmse_calibration = sqrt(1./(size(rws,1)-1).*sum((water_compare-water).^2))
			rmse_ratio_ws = rmse_calibration/100*(rmse_calibration/(100-rmse_calibration)+1)

			figure()
			scatter(rws[:,3],water[:]./(100.-water[:]))
			plot(rws_calibration,water_ratio_calibration,color="red",linewidth=2.0)
			plot(rws_calibration,water_ratio_calibration+rmse_ratio_ws,color="red",linestyle="--",linewidth=1.0)
			plot(rws_calibration,water_ratio_calibration-rmse_ratio_ws,color="red",linestyle="--",linewidth=1.0)
			plot(rws_calibration,water_ratio_calibration,color="red")
			xlim(0,maximum(rws_calibration)+1./4*maximum(rws_calibration))
			ylim(0,maximum(water_ratio_calibration)+1./4*maximum(water_ratio_calibration))

			xlabel(L"A$_{water}$/A$_{silicates}$, area ratio",fontsize=18,fontname="Arial")
			ylabel("Water/Glass, weight ratio",fontsize=18,fontname="Arial")
			annotate(L"R$_{ws}$=",xy=(0.3,0.9),xycoords="axes fraction",fontsize=18,fontname="Arial",horizontalalignment="center")
			annotate("$(round(coef[1],5)) +/- $(round(sigma[1],5))",xy=(0.3,0.8),xycoords="axes fraction",fontsize=18,fontname="Arial",horizontalalignment="center")
			annotate("Standard deviation\n= $(round(rmse_calibration,2)) wt%",xy=(0.7,0.3),xycoords="axes fraction",horizontalalignment="center",fontsize=18,fontname="Arial")
			savefig(paths[6])

			writecsv(string(paths[5]), ["spectrum" "product" "Rws" "Water input" "Water Raman";liste[:,1] liste[:,2] rws[:,3] water water_compare])
			#return [rws[:,3] water water_compare]
		else
			water_predicted = 100.*(rws[:,3].*prediction_coef[1]./(rws[:,3].*prediction_coef[1]+1))
			water_predicted_high = 100.*(rws[:,3].*sum(prediction_coef)./(rws[:,3].*sum(prediction_coef)+1))
			writecsv(string(paths[5]), ["spectrum" "product" "Rws" "Water Raman" "Error";liste[:,1] liste[:,2] rws[:,3] water_predicted water_predicted_high-water_predicted])
			#return [rws[:,3] water_predicted]
		end

	elseif switches[1] == "external"
		liste=readcsv(paths[1], skipstart=1,use_mmap = mmap_switch)

		reference = liste[:,1]
		# reference product names are in first column for record
		reference_density = liste[:,3]
		reference_water = liste[:,4]

		sample = liste[:,5]
		# sample product names are in column 6 for record
		sample_density = liste[:,7]

		water = ones(size(liste,1),1) # the array to record water content

		for i = 1:size(liste,1)

			#importing data
			reference_sp = readdlm(string(paths[2],reference[i]),input_properties[1],skipstart=input_properties[2],use_mmap = mmap_switch)
			sample_sp = readdlm(string(paths[2],sample[i]),input_properties[1],skipstart=input_properties[2],use_mmap = mmap_switch)

			if reference_sp[end,1] < reference_sp[1,1]
				reference_sp = flipdim(reference_sp,1)
			end

			if sample_sp[end,1] < sample_sp[1,1]
				sample_sp = flipdim(sample_sp,1)
			end

			# temperature - excitation line correction with the Hehlen formula
			x_ref, reference_sp_long, ~ = tlcorrection(reference_sp,temperature,laser,correction="hehlen",normalisation="no",density=reference_density[i])
			x_sp, sample_sp_long, ~ = tlcorrection(sample_sp,temperature,laser,correction="hehlen",normalisation="no",density=sample_density[i])

			# Linear baseline subtraction
			reference_sp_corr2, reference_baseline = baseline(x_ref,reference_sp_long,roi_hf_external,"poly",p=1.0)
			sample_sp_corr2, sample_baseline = baseline(x_sp,sample_sp_long,roi_hf_external,"poly",p=1.0)

			# Area calculation between 3100 and 3800 cm-1
			Area_reference = trapz(x_ref[roi_hf_external[1,2] .<x_ref[:,1].<roi_hf_external[2,1]],reference_sp_corr2[roi_hf_external[1,2] .<x_ref[:,1].<roi_hf_external[2,1]])
			Area_sample = trapz(x_sp[roi_hf_external[1,2] .<x_sp[:,1].<roi_hf_external[2,1]],sample_sp_corr2[roi_hf_external[1,2] .<x_sp[:,1].<roi_hf_external[2,1]])

			# With the intensity
			Intensity_reference = maximum(reference_sp_corr2[roi_hf_external[1,2] .<x_ref[:,1].<roi_hf_external[2,1]])
			Intensity_sample = maximum(sample_sp_corr2[roi_hf_external[1,2] .<x_ref[:,1].<roi_hf_external[2,1]])

			# Water calculation
			water[i] = (reference_water[i]./1.8*reference_density[i]) .*Area_sample./Area_reference # water is in mol/L
			#water[i] = (reference_water[i]./1.8*reference_density[i]) .*Intensity_sample./Intensity_reference # water is in mol/L

			water[i] = water[i] *1.8./sample_density[i] # water converted in wt%

			figure(figsize=(7.5,15))

			subplot(311)
			title("External calibration: $(round(water[i],2)) wt% water,\n sp. $(sample[i]) with ref. $(reference[i]) ")
			plot(reference_sp[:,1],reference_sp[:,2],color="black",label="Reference")
			plot(sample_sp[:,1],sample_sp[:,2],color="blue",label="Sample")
			xlabel(L"Raman shift, cm$^{-1}$",fontsize=18,fontname="Arial")
			legend(loc="best",frameon=false)

			subplot(312)
			plot(x_ref[:],reference_sp_long[:],color="black",label=L"T-$\nu$ corrected reference")
			plot(x_ref[:],reference_baseline,color="grey",label="Reference baseline")
			plot(x_sp[:],sample_sp_long[:],color="blue",label=L"T-$\nu$ corrected sample")
			plot(x_sp[:],sample_baseline,color="cyan",label="Sample baseline")
			xlabel(L"Raman shift, cm$^{-1}$",fontsize=18,fontname="Arial")
			ylabel("Normalized intensity, a. u.",fontsize=18,fontname="Arial")
			legend(loc="best",frameon=false)

			subplot(313)
			plot(x_ref,reference_sp_corr2,color="black",label="Ref. spectrum")
			plot(x_sp,sample_sp_corr2,color="blue",label="Sample spectrum")
			xlabel(L"Raman shift, cm$^{-1}$",fontsize=18,fontname="Arial")
			legend(loc="best",frameon=false)
			tight_layout()
			# saving the water content
			writecsv(string(paths[5]), ["Spectrum" "Sample Name" "Estimated water content";liste[:,5] liste[:,6] water])

			#Saving the generated figure (specified directory should exist)
			savefig(string(paths[4],"External_Std_",sample[i],".pdf"))

		end


	else # or error message
		error("The first switch indicates if you want to use an internal or external calibration mode. Please choose between both.")
	end
end
