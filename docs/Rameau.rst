.. _Rameau:
***********************
RamEau
***********************

-------------------
Introductory notes
-------------------

This is the Julia version of the RamEau software. It allows quantification of the water content of glasses following the protocol described in 

Le Losq, Neuville, Moretti, Roux, 2012. Determination of water content in silicate glasses using Raman spectrometry: Implications for the study of explosive volcanism. American Mineralogist 97, 779-790.

The Rameau Pascal/fortran initial software is available through the american mineralogist website.

The function rameau allows you to calculate the area ratio between the water and silicate signals, and to relate it to measured water content for establishing a new calibration, or to predict the glass water content with the knowledge of the calibration factor.

This version of Rameau aims to bring new features with a lot of flexibility over the way you are calculating the calibration.

The current version mostly work following the protocol described in Le Losq et al. (2012). However, I hope bringing new features in the future such as a double baseline feature:

- A first baseline is subtracted from the raw data under the water peak, giving a first corrected spectrum;
- A second one is subtracted from the first corrected spectrum.

This is purely experimental for now so I don't advise using it for now, but we have reasons to perform such approach, such as:

- It may be better (or not...) to remove a general spectral background before performing the Long correction;
- Different smoothing coefficients for the silicate and water region may provide more robust and precise baseline fitting.

Please read carefully the following description, and after that jump into the examples section of Spectra.jl to see Spectra.rameau in action on the dataset published in 2012! To conclude, any bug report, contributions on Github and suggestions will help improving this software and Spectra.jl in general. So you're very welcome to provide any feedback!

NOTE ON ABREVIATIONS: Rws in the following refers to the ratio between the area of the water peak and that of the silicate bands.

------------------------------
Function rameau
------------------------------

Call rameau as

    rameau(paths::Tuple,input_properties::Tuple,switches::Tuple;prediction_coef=0.035,temperature=23.0,laser=532.0,lb_break=2010.0,hb_start=1000.0)

INPUTS:
	
    paths: Tuple{Strings}, it contains the following strings: 
	
		in_liste: the relative path of your liste of spectra (see the description of this liste below), e.g. "./liste_2012.csv"
		
		in_path = the relative path of the repertory where your raw spectra are stored, e.g. "./raw/"
		
		out_path = the relative path of the repertory where you want to output the corrected spectra, e.g. "./treated/" (THIS FOLDER MUST EXIST PRIOR TO RUN)
		
		fig_path= the relative path of the repertory for saving the figures, e.g. "./figures/"
		
		rws_save_file = the relative path of the csv file where the results will be stored, e.g. "rws_calculated.csv". If calibration is set to "yes", it outputs an array containing the Rws, the provided water content and the calculated water content in columns 1, 2 and 3, respectively. If calibration is set to anything else, it outputs the rws and the water content predicted with the provided prediction_coef
		
		rws_save_fig = the relative path of the figure for the calibration, e.g. "rws_calculated.jpg". Note: It only works when calinration is on "yes"
		
	input_properties: Tuple, it contains the delimiter of the columns in your spectra files, and the number of starting lines to skip. For instance, use ('\t',1) if your files have one header line and the columns are separated  by tabulations, or (',',0) if you use CSV files without header.
	
	switches: Tuple, it contains the different switches as (type of calibration, calibration?, double baseline, temperature-laser wavelength correction). 
		
		Type of calibration: should be set to "internal" or "external". The external mode requires standards, as described in Behrens et al. (2006) and Thomas et al. (2008). If you use the external mode, you can leave blanks for the other switches. For instance, enter ("external","","","").
		
		calibration? :  use only with the "internal" mode. Enter yes if you are setting up a new calibration, or "no" if you want to use a previously determined calibration coefficient with providing the later in the prediction_coef variable (see options).
		
		double baseline? : use only with the "internal" mode. Enter "yes" if you want to use it. This is a new feature for internal calibration, see comments below.
		
		temperature-laser wavelength correction? : use only with the "internal" mode. Enter "yes" if you want to use it. This asks if you want to use the temperature-laser wavelength correction as done in Le Losq et al. (2012). If you use the double baseline mode, this correction is applied after removing the background under the water peak.
	
	For instance, you have a bunch of spectra and you want to perform a new calibration, without the double baseline but the Long correction (you are following the "traditional approach published in 2012). So you are going to enter ("internal",""yes","no","yes"). Now, you think removing the initial background may be good before calculating the calibration. So you will enter ("internal","yes","yes","yes"). Now that you have determined your calibration coefficient, you can make any prediction on a bunch of new spectra in the futur by entering ("internal","no","yes","yes") and by providing the good calibration coefficient in the variable prediction_coef (see options below).
	
OPTIONS:
	
	prediction_coef: Float64, this is the calibration coefficient that will be used in the predictions, if you use the predictive mode (i.g., calibration switch is NOT set to "yes"). Default = 0.069.
	
	temperature: Float64, the temperature for the temperature-laser wavelength correction in Celsius. Default = 23.0.
	
	laser: Float64, the laser wavelength in nm for the temperature-laser wavelength correction. Default = 532.0.
	
	lb_break: Float64, for double baseline correction, the breaking point before which the software will consider the BIRs in the low frequency region. Default = 2010.0.
	
	hb_start: Float64, for double baseline correction, the breaking point after which the software will consider the BIRs in the high frequency region. Default = 1000.0.
	
	basetype: String, the type of baseline you want to fit. Corresponds to the "basetype" parameter of the baseline function.  Default = "gcvspline".
	
OUTPUTS:

	For now reameau does not provide any outputs, but save everything in the files you indicate in the variable "paths". This may change in the futur, so stay tuned!

-----------------------------------
Note on the input file liste
-----------------------------------

The great news about RamEau in Julia is that you can work your file liste in Excel, as it is now a CSV file. It makes it much more robust and readable.

If using the "internal" mode, this file liste MUST contain:

	column 1: the file name and extensions, e.g. myspectrum.txt

	column 2: the name of your product

	column 3: the water content, if known. If unknow, put 0.0

	column 4: the spline coefficient for the silicate part. Note: this value is used in the single baseline procedure for the whole spectrum

	column 5: the spline coefficient for the water part, in case you use the double baseline fitting procedure

	columns 6 to end: the beginning and ends of the BIRs, paired. Please keep the same number of BIRs for all the spectra in one batch.

If using the "external" mode, this file liste MUST contain:

	column 1: the file name and extensions of the references, e.g. myreference.txt

	column 2: the name of your references

	column 3: the water content of the references, in wt%

	column 4: the density of the references, in kg m-3
	
	column 5: the file name and extensions of the samples, e.g. mysample.txt

	column 6: the name of your samples

	column 7: the estimated density of your samples, in kg m-3

WARNING: BE SURE THAT THE NUMBER YOU PROVIDE ARE FLOAT NUMBER!

-----------------------------------------------------------------------------
Note on the temperature and excitation line effects corrections
-----------------------------------------------------------------------------

The "internal" mode uses the "long" mode of the tlcorrection function, whereas the "external" mode uses the "hehlen", which takes into account the sample density (see tlcorrection function documentation). This allows to intrisically correct the intensity from density effects.

-----------------------------------
Note on the double baseline feature
-----------------------------------

This is purely experimental and will probably strongly change in the upcoming future. However, some thoughts about why we may enjoy such function:

This is to be used with the internal calibration mode.

I added this step to avoid the strong distortion of the spectra during the Long correction. Indeed, spectra are distorded because even the parts without signals are not close to a zero intensity in raw spectra. Therefore, to avoid that, I added this feature which basically fits a linear function between 1300 and 2000 cm-1, where no signals are usually expected in silicate glasses. I futher take the shot to fit the water peak at the same time. Then, a second baseline will fit the basis of the peaks below 1500 cm-1. This double baseline approach allows to avoid a strong distortion of the signal due to the Long correction, and further allow working with different spline coefficients.

However, I warn the user that this is not always the best solution... Indeed, the slight signal distortion created by the usual Long correction sometimes helps fitting the baseline, as it nearly create a flat, linear increase of the spectral background.

Therefore, this is up to the user to choose what is best in his case.

From my test, switching from one mode to the other might improve or worsen the standard deviation of the calibration of around 0.1-0.3 wt%. It might (or not...) improve the robustness of the baseline fitting procedure.

-----------------------------------------------------------------------
Note on the use of KRregression baseline fitting instead of GCV splines
-----------------------------------------------------------------------

This is to be used with the internal calibration mode.

Back in 2012 we mostly used the Generalized Cross-Validated splines for fitting the spectral background. However, recent developments show that KRregression or SVMregression may provid better results with less headache for the user (not need to tune the spline coefficient parameter). From experience, using a spline carefully adjusted provides better result. However, using KRregression may provide good results without headache to adjust any parameter. For now this is an experimental feature.