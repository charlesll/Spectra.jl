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
# This files contains functions to remove the signal of a cristal from the glass Raman signal.
#
#
#############################################################################

"""
	ctxremoval(liste,in_path,out_path,roi_all;input_properties=('\t',0),algorithm="FastICA",plot_intermediate_show = "no",plot_mixing_show = "yes",plot_final_show = "no",save_fig_switch = "yes", shutdown = 1300.,scaling=100.)
	
INPUTS

	liste: Array{Float64}, an array contaning the information for getting the spectra. 
	
		Column 1: the name and relative path for the crystal spectra; 
		
		Column 2: the name and relative path for the (mixed) glass spectra;
		
		Column 3: smo_water, the coefficient of smothing for the spline that fits the backgrous below the signal of water (if any) in the glass spectra
		
		Column 4: number of iteration for the generation of new mixed spectra
		
		Column 5: the K_Start parameter, set at 0.0
		
		Column 6: the K_increament parameter, for mixing the signals

	in_path: String, the relative location of the data, e.g. "./raw/"

	out_path: String, the relative location where you want to save the corrected spectra, e.g. "./treated/" 
	
	roi_all: Tuple, contains 2 arrays and 2 Float64 numbers. The 2 arrays indicate the regions of interest where the background correction is applied, for the cristal and the glass. The  2 float numbers indicate the starting and ending frequency of the peak used to correct the spectra from any shift in frequency. For instance:
	
		Those are the roi for fitting the baseline on the crystal (roi_ctx) and glass (roi_glass) signals:
		
		roi_ctx = [1260. 2000.;2000. 4000.]
		
		roi_glass = [1260. 2000.;2000. 3000.;3750. 4000.]
		
		We have a strong peak from the crystal at ~650 cm-1 that we can use to correct the spectra from any shift in frequency. So we indicate the values here:
		
		roi_xshift_low = 655.
		
		roi_xshift_high = 670.
		
		Then we construct the final tuple as:
		
		roi_all = (roi_ctx,roi_glass,roi_xshift_low,roi_xshift_high)

OPTIONS

	input_properties: Tuple, this tuple contains the delimiter and the number of lines to skip in the raw data files. Default = ('\t',0);
	
	algorithm: String, This indicates if the FastICA or the Non-negative Matrix Factorisation (NMF) algorithms from SciKit Learn will be used. Default = "FastICA";
	
	plot_intermediate_show: String, This should be equal to "yes" or "no". It displays the intermediate figures. Default = "no";
	
	plot_mixing_show: String, This should be equal to "yes" or "no". It displays the figures showing the mixing step. Default = "yes";
	
	plot_final_show: String, This should be equal to "yes" or "no". It displays the final figures, showing the background subtraction and the retrieved signals. Default = "yes";

	save_fig_switch: String, This should be equal to "yes" or "no". It indicates if you want to save the final figures in the location indicated by out_path;
	
	shutdown: Float64, indicates where you consider the signals from silicate units to stop. Default = 1300.0;
	
	scaling: Float64, the retrieved spectra are scaled to the original spectra using the Boson peak, located ~ 60-80 cm-1. This parameters indicates where you consider the Boson peak to stop for the scaling procedure. No need to put a too high value, as you might get strong crystal signals at frequencies > 100-150 cm-1.

OUTPUTS

	All the corrected spectra and figures are saved in the location indicated in out_path. No direct outputs in Julia.
"""
function ctxremoval(liste,in_path,out_path,roi_all;input_properties=('\t',0),algorithm="FastICA",plot_intermediate_show = "no",plot_mixing_show = "yes",plot_final_show = "no",save_fig_switch = "yes", shutdown = 1300.,scaling=100.)
	
	# Model for final adjustement
	model_ctx_adj(x, p) = p[1]*x
	
	# GRABING THE REGIONS OF INTEREST
	roi_ctx = roi_all[1]
	roi_glass = roi_all[2]
	roi_xshift_low = roi_all[3]
	roi_xshift_high = roi_all[4]
	
	# NUMBER OF FILES?
	nb_exp = size(liste,1)

	# SARTING THE LOOP TO TREAT EACH COUPLE OF SPECTRA
	for i = 1:nb_exp
    
	    # IMPORTING DATA
	    ctx = readdlm(in_path*liste[i,1],input_properties[1],skipstart=input_properties[2])
	    glass = readdlm(in_path*liste[i,2],input_properties[1],skipstart=input_properties[2])
    
	    #CHECKING THAT X IS INCREASING
	    if glass[end,1] < glass[1,1]
	        glass = flipdim(glass,1)
	    end
    
	    if ctx[end,1] < ctx[1,1]
	        ctx = flipdim(ctx,1)
	    end
    
	    # THE COMMON X AXIS, TAKEN AS THE GLASS ONE
	    true_x = collect(minimum(glass[:,1]):0.8: maximum(glass[:,1]))
        
	    # RESAMPLING THE TWO SPECTRA WITH THE SAME X
	    data = zeros(size(true_x,1),2)
    
	    spl_ctx = Spline1D(ctx[:,1],ctx[:,2])
	    spl_gls = Spline1D(glass[:,1],glass[:,2])
    
	    data[:,1] = evaluate(spl_ctx, true_x)
	    data[:,2] = evaluate(spl_gls, true_x)
    
	    # A VARIABLE TO STORE THE BASELINE CORRECTED SPECTRA
	    bas_matrix = zeros(size(true_x,1),2)
    
	    # WE GET THE SMOOTHING SPLINE COEFFICIENT
	    smo_water = liste[i,3]
    
	    # BASELINE CTX
	    ycorr_ctx_1, baseline1_ctx = baseline(true_x[true_x[:].<roi_ctx[1,2]], data[true_x[:].<roi_ctx[1,2],1],reshape(roi_ctx[1,:],1,2),"poly",p=1.0)
	    ycorr_ctx_2, baseline2_ctx = baseline(true_x[true_x[:].>=roi_ctx[1,2]],data[true_x[:].>=roi_ctx[1,2],1],reshape(roi_ctx[2,:],1,2),"gcvspline",p=smo_water)
	    baseline_ctx = [baseline1_ctx;baseline2_ctx]
	    bas_matrix[:,1] = [ycorr_ctx_1;ycorr_ctx_2]
    
	    # BASELINE GLASS
	    ycorr_gls_1, baseline1_gls = baseline(true_x[true_x[:].<roi_glass[1,2]], data[true_x[:].<roi_glass[1,2],2],reshape(roi_glass[1,:],1,2),"poly",p=1.0)
	    ycorr_gls_2, baseline2_gls = baseline(true_x[true_x[:].>=roi_glass[1,2]],data[true_x[:].>=roi_glass[1,2],2],roi_glass[2:end,:],"gcvspline",p=smo_water)
	    baseline_gls = [baseline1_gls;baseline2_gls]
	    bas_matrix[:,2] = [ycorr_gls_1;ycorr_gls_2]
    
	    # DISPLAYING INTERMEDIATE PLOTS
	    if plot_intermediate_show == "yes"
	        figure()
			
	        title("Raw spectra and baseline")
	        xlabel(L"Raman shift, cm$^{-1}$")
	        ylabel("Intensity, counts")

	        plot(true_x,data[:,1],color="black",label="Ctx file "liste[i,1])
	        plot(true_x,data[:,2],color="blue",label="Glass file "liste[i,2])
        
	        plot(true_x,baseline_ctx,color="red")
	        plot(true_x,baseline_gls,color="orange")
	        legend(loc="best",fancybox=true)
			show()
			
	        figure()
			
	        title("Raw spectra and baseline")
	        xlabel(L"Raman shift, cm$^{-1}$")
	        ylabel("Intensity, counts")

	        plot(true_x,data[:,1],color="black",label="Ctx file "liste[i,1])
	        plot(true_x,data[:,2],color="blue",label="Glass file "liste[i,2])
        
	        plot(true_x,baseline_ctx,color="red")
	        plot(true_x,baseline_gls,color="orange")
	        legend(loc="best",fancybox=true)
			xlim(2600,4000)
			ylim(0,maximum(data[2600 .<true_x.<4000,:]))
	        show()
	    end
    
	    if plot_intermediate_show == "yes"
	        figure()
	        title("Frequency shift correction")
	        			plot(true_x[roi_xshift_low.<true_x.<roi_xshift_high],bas_matrix[roi_xshift_low.<true_x.<roi_xshift_high,1]./maximum(bas_matrix[roi_xshift_low.<true_x.<roi_xshift_high,1]),color="blue",label="shifted")
	        			plot(true_x[roi_xshift_low.<true_x.<roi_xshift_high],bas_matrix[roi_xshift_low.<true_x.<roi_xshift_high,2]./maximum(bas_matrix[roi_xshift_low.<true_x.<roi_xshift_high,2]),color="red",label="reference")
	    end
    
	    # AUTOMATIC ADJUSTEMENT OF ANY X SHIFT BETWEEN THE SPECTRA, the glass spectrum is the reference
	    ~, bas_matrix[:,1], p = xshift_correction(true_x,bas_matrix[:,1],true_x[roi_xshift_low .< true_x .<roi_xshift_high], data[roi_xshift_low .< true_x .<roi_xshift_high,2],data[roi_xshift_low .< true_x .<roi_xshift_high,1])
    
	    # DISPLAYING INTERMEDIATE PLOTS
	    if plot_intermediate_show == "yes"
	        			plot(true_x[roi_xshift_low.<true_x.<roi_xshift_high],bas_matrix[roi_xshift_low.<true_x.<roi_xshift_high,1]./maximum(bas_matrix[roi_xshift_low.<true_x.<roi_xshift_high,1]),color="green",label="corrected")
	        			plot(true_x[roi_xshift_low.<true_x.<roi_xshift_high],bas_matrix[roi_xshift_low.<true_x.<roi_xshift_high,2]./maximum(bas_matrix[roi_xshift_low.<true_x.<roi_xshift_high,2]),color="red",label="reference")
	        legend(loc="best",fancybox=true)
	        xlabel(L"Raman shift, cm$^{-1}$")
	        ylabel("Intensity, counts")
	    	show()
			
	        figure()
	        title("Spectra after the baseline and X shift corrections")
	        plot(true_x,bas_matrix[:,1],color="black",label="Cristal file "liste[i,1])
	        plot(true_x,bas_matrix[:,2],color="blue",label="Glass file "liste[i,2])
	        xlabel(L"Raman shift, cm$^{-1}$")
	        ylabel("Intensity, counts")
	        legend(loc="best",fancybox=true)
			show()
	    end
    
	    ######################################################################
	    # CRYSTAL SIGNAL REMOVAL ALGORITHM
	    # The following code creates new spectra with mixed glass-ctx signals
	    # Those spectra are used in FastICA to recovert the pure glass and ctx components
    
	    # READING THE GOOD VALUES FOR THE ITERATIVE PROCESS
	    nb_iter = liste[i,4] # the number of iterations
	    K_start = liste[i,5] # for the iteration, starting mixing point (fraction of ctx signal)
	    K_increment = liste[i,6] # for the iteration, increamentation mixing point (fraction of ctx signal)
    
	    # CREATING USEFULL FUNCTIONS
	    function f1(k,ctx,g_c)
	        g = (g_c - k.*ctx)./(1.0-k)
	        return g
	    end
    
	    function f2(k,ctx,g) # maybe we do not need this one ?
	        g_c = k.*ctx+g
	        return g_c
	    end
    
	    # WE ARE ONLY RECODING SIGNALS BELOW 1300 cm-1 (strictly)
	    x_for_ica = true_x[true_x.<shutdown]
	    y_for_ica = ones(size(x_for_ica,1),nb_iter)
    
	    # LAST FEW THINGS BEFORE ITERATIONS
	    results = zeros(size(bas_matrix)) # the final matrix with mixed signals
	    K = K_start # we start with K_start ctx signal
    
	    for iter = 1:nb_iter # iterative removal, should be good in a few iterations
	        if iter == 1
	            results[:,:] = bas_matrix[:,:]
	            results[:,1] = K*bas_matrix[:,1]
	        else
	            results[:,1] = K.*bas_matrix[:,1]
	            results[:,2] = f1(K,bas_matrix[:,1],bas_matrix[:,2])
				#results[:,2] = f2(K,bas_matrix[:,1],bas_matrix[:,2])
	            results[true_x.>shutdown,2] = bas_matrix[true_x.>shutdown,2]
	        end
        
	        if plot_mixing_show == "yes" # for showing intermediate plots
				if iter == 1
				    figure()
				    title("The generation of spectra for FastICA")
		            subplot(311)
		            title("Ctx, cristal signals",fontsize=18,fontname="Arial")
		            plot(true_x,results[:,1],color=[iter/(nb_iter+0.0),0.,0.])
		            xlim(0,shutdown)
		            subplot(312)
		            title("G, unmixed glass signals = (G_C - K*Ctx)/(1-K)",fontsize=18,fontname="Arial")
		            ylabel("Intensity, counts",fontsize=18,fontname="Arial")
		            xlim(0,shutdown)
		            plot(true_x,results[:,2],color=[iter/(nb_iter+0.0),0.,0.])
		            subplot(313)
		            title("G_C, mixed glass signals = G + K*CTX",fontsize=18,fontname="Arial")
		            plot(true_x,results[:,1]+(1-K).*results[:,2],color=[iter/(nb_iter+0.0),0.,0.])
		            xlabel(L"Raman shift, cm$^{-1}$",fontsize=18,fontname="Arial")
		            xlim(0,shutdown)
					tight_layout()
				else
					subplot(311)
					plot(true_x,results[:,1],color=[iter/(nb_iter+0.0),0.,0.])
					subplot(312)
					plot(true_x,results[:,2],color=[iter/(nb_iter+0.0),0.,0.])
					subplot(313)	
					plot(true_x,results[:,1]+(1-K).*results[:,2],color=[iter/(nb_iter+0.0),0.,0.])
				end
	            
	        end
        
			y_for_ica[:,iter] = results[true_x.<shutdown,2] # old way: actually works well
			#y_for_ica[:,iter] = K*bas_matrix[true_x.<shutdown,1] + bas_matrix[true_x.<shutdown,2]# New trial: only positive addition: may improve NMF
		
	        # INCREMENTING K
	        K = K + K_increment
        
			
	    end
		if plot_mixing_show == "yes"
			show()
		end
    
		# We test the algorithm wanted by the user
		if algorithm == "FastICA"
			
		    # PREPROCESSING (STANDARDIZATION)
		    X_scaler = preprocessing[:StandardScaler]()
		    X_scaler[:fit](y_for_ica)
		    y_for_ica_sc = X_scaler[:transform](y_for_ica)
			
	    	# USING FastICA ON THE G SIGNALS
	    	model = decomposition[:FastICA](n_components=2, fun = "cube")
	    	S = model[:fit_transform](y_for_ica_sc)
    
	    	# TREATING THE OUTPUTS, FastICA can return negative components so we treat the S array for that
	    	S_corr = S[:,:]
	    	for col = 1:2
	        	if maximum(S[990. .< x_for_ica .< 1050,col]) < maximum(S[1200 .< x_for_ica .<shutdown,col]) # checking if something is negative
	            	S_corr[:,col] = -S[:,col]
	        	end
	        	S_corr[:,col] = (S_corr[:,col]-minimum(S_corr[:,col])) / maximum(S_corr[:,col]-minimum(S_corr[:,col])) # normalisation
	    	end
		elseif algorithm == "NMF"
			
		    # PREPROCESSING (STANDARDIZATION)
		    X_scaler = preprocessing[:MinMaxScaler]()
		    X_scaler[:fit](y_for_ica)
		    y_for_ica_sc = X_scaler[:transform](y_for_ica)
			
			# with sklearn
	    	#model = decomposition[:NMF](n_components=2,init="nndsvda")
	    	#S_corr = model[:fit_transform](y_for_ica_sc)
			
			# with julia NMF implementation: nndsvda init + alspgrad algo
			r = nnmf(transpose(y_for_ica), 2; init=:nndsvda, alg=:alspgrad, maxiter=100, verbose=false)			 
			H = r.H; W = r.W; S_corr = transpose(H./maximum(H,2))
			
		else
			error("Not implemented, choose between NMF and FastICA")
		end
		
	    # WHICH SIGNAL IS WAHT? Now we determine who is the glass? => the Boson peak is more intense!
	    if maximum(S_corr[40 .< x_for_ica .< 100,2]) < maximum(S_corr[40 .< x_for_ica .<100,1])
	        glass_ica = S_corr[:,1] - minimum(S_corr[:,1])
	        ctx_ica = S_corr[:,2] - minimum(S_corr[:,2])
	    else 
	        glass_ica = S_corr[:,2] - minimum(S_corr[:,2])
	        ctx_ica = S_corr[:,1] - minimum(S_corr[:,1])
	    end
    	
	    # FINAL RECORDING MATRIX
		fit = curve_fit(model_ctx_adj, glass_ica[x_for_ica.<scaling], bas_matrix[true_x .<scaling,2], [1.0])
		coef_scale = fit.param
	    spectra_final = [true_x [model_ctx_adj(glass_ica,coef_scale);bas_matrix[true_x.>=shutdown,2]]]
		
	    # FINAL FIGURE, SHOWED
	    figure(figsize=(15,15))
		
	    subplot(3,2,(1,2)) # UPPER PLOT
		title("Glass sp. no $(liste[i,2]), Baseline correction",fontsize=18,fontname="Arial")
	    plot(true_x,bas_matrix[:,1],color="blue",label="Cristal")
	    plot(true_x,bas_matrix[:,2],color="black",label="Glass")
	    ylabel("Intensity, counts",fontsize=18,fontname="Arial")
	    xlim(0,4000)
	    legend(loc="best",fancybox=true)
		
	    subplot(3,2,3)
	    title("results from the algorithm "*algorithm,fontsize=18,fontname="Arial")
	    plot(x_for_ica,glass_ica,color="black")
	    plot(x_for_ica,ctx_ica,color="blue")
	    ylabel("Intensity, counts",fontsize=18,fontname="Arial")
	    xlim(0,shutdown)
    
	    subplot(3,2,4)
	    title("Final fit, silicate region",fontsize=18,fontname="Arial")
	    plot(x_for_ica, model_ctx_adj(glass_ica,coef_scale),color="green",label="Glass, corrected")
	    plot(x_for_ica, bas_matrix[true_x.<shutdown,2],color="black",label="Glass, initial")
	    xlabel(L"Raman shift, cm$^{-1}$",fontsize=18,fontname="Arial")
	    ylabel("Intensity, counts",fontsize=18,fontname="Arial")
	    xlim(0,shutdown)
    
	    subplot(3,2,(5,6)) # LOWER PLOT
		title("Final fit, entire spectra",fontsize=18,fontname="Arial")
		plot(true_x,bas_matrix[:,2],color="black",label="Glass, initial")
	    plot(spectra_final[:,1],spectra_final[:,2],color="green",label="Glass, corrected")
	    legend(loc="best",fancybox=true)
	    xlabel(L"Raman shift, cm$^{-1}$",fontsize=18,fontname="Arial")
	    ylabel("Intensity, counts",fontsize=18,fontname="Arial")
	    xlim(0,4000)
		
		tight_layout()
	    # AND IF THE USER WANTS TO SAVE THE FIGURE, WE OUTPUT THEM ALSO IN OUT_PATH
	    if save_fig_switch == "yes"
	        savefig(out_path*liste[i,2]*".pdf")
	    end
		
		if plot_final_show == "yes"
			show()
		else
			close()
		end
	
	    # WE WRITE THE CORRECTED SPECTRA IN THE CORRECT OUTPUT PATH (OUT_PATH)
	    writedlm(out_path*liste[i,2],spectra_final)
    
	end
end
	
	