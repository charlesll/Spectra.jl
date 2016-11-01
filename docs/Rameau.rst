.. _Rameau:
***********************
RamEau
***********************

-------------------
Introductory notes
-------------------

This is the Julia version of the RamEau software. It allows quantification of the water content of glasses following the internal and external protocols described in:

Thomas, R. 2000. “Determination of Water Contents of Granite Melt Inclusions by Confocal Laser Raman Microprobe Spectroscopy.” American Mineralogist 85 (5-6): 868–72.

Behrens, Harald, Jacques Roux, Daniel R. Neuville, and Michael Siemann. 2006. “Quantification of Dissolved H2O in Silicate Glasses Using Confocal microRaman Spectroscopy.” Chemical Geology 229 (1-3): 96–112. doi:10.1016/j.chemgeo.2006.01.014.

Le Losq, Neuville, Moretti, Roux, 2012. Determination of water content in silicate glasses using Raman spectrometry: Implications for the study of explosive volcanism. American Mineralogist 97, 779-790.

The Rameau Pascal/fortran initial software is available through the american mineralogist website. This version goes much beyond the previous version. It allows using various modes for internal calibration, and further allows using external calibrations too.

Internal calibration mode refers to the technic of using the silicate peaks to scale the water peak, before relating this ratio to the sample water concentration. External calibrations directly refer the integrated intensity of peak height of the O-H stretching band to the water content, through the use of a standard glass for which this relationship is well constrained. It assumes a linear relationship between the water peak height and the glass water content. See the references listed above for more details.

Please read carefully the following description, and after that jump into the examples section of Spectra.jl to see Spectra.rameau in action on a fraction of the dataset published in 2012. For the full dataset, please consult the American Mineralogist website. To conclude, any bug report, contributions on Github and suggestions will help improving this software and Spectra.jl in general. So you're very welcome to provide any feedback!

NOTE ON ABREVIATIONS: Rws in the following refers to the ratio between the area of the water peak and that of the silicate bands.

------------------------------
Function rameau
------------------------------

Call rameau as

    rameau(paths::Tuple,input_properties::Tuple,switches::Tuple;prediction_coef=0.035,temperature=23.0,laser=532.0,lb_break=2010.0,hb_start=1000.0,mmap_switch=true)

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

	prediction_coef: Float64, this is the calibration coefficient that will be used in the predictions, if you use the predictive mode (i.g., calibration switch is NOT set to "yes"). Default = 0.069.

	temperature: Float64, the temperature for the temperature-laser wavelength correction in Celsius. Default = 23.0.

	laser: Float64, the laser wavelength in nm for the temperature-laser wavelength correction. Default = 532.0.

	lb_break: Float64, for double baseline correction, the breaking point before which the software will consider the BIRs in the low frequency region. Default = 2010.0.

	hb_start: Float64, for double baseline correction, the breaking point after which the software will consider the BIRs in the high frequency region. Default = 1000.0.

	basetype: String, the type of baseline you want to fit. Corresponds to the "basetype" parameter of the baseline function.  Default = "gcvspline".

  mmap_switch: false or true, this allows to switch on or off the memory mapping in the `readcsv`/`readdlm` functions that `rameau` uses. Default = "true".

OUTPUTS:

	Rameau does not provide any outputs directly in Julia, but saves everything in the folders you indicate in the variable "paths".

--------------
Quick examples
--------------

In this example, the Julia code and the csv liste (myliste.csv) of spectra are in the working folder, the data are in ./raw/, and we want to output the corrected spectra and the figures in the ./treated/ and ./figures/ folders. So we set things like:

	in_liste: "./myliste.csv"

	in_path = "./raw/"

	out_path = "./treated/"

	fig_path= "./figures/"

	rws_save_file = "./treated/"

	rws_save_fig = "./figures/mycalibration.pdf"

	paths = (in_liste,in_path,out_path,fig_path,rws_save_file,rws_save_fig)

Now, for performing an internal calibration as explained in Le Losq et al. (2012), enter:

	switches = ("internal",""yes","no","yes")

and call Rameau:

	rameau(paths,switches,input_properties = ('\t',0))

This will allow you to get your prediction coefficient prediction_coef With this knowledge, you can predict values from the spectra of new glasses with the names in "myliste_newglasses.csv" with using the commands:

	in_liste = "myliste_newglasses.csv"

	switches = ("internal",""no","no","yes")

	rameau(paths,switches,prediction_coef = 0.0059, input_properties = ('\t',0))

For an external calibration, you need a standard glass with known water concentration. You also need the knowledge of the densities of the standard and sample glasses. Then, the following commands allow you to calculate the water content of your sample with using the protocol described in Thomas et al. (2008; see also references cited therein):

	in_liste: "./myliste.csv"

	in_path = "./raw/"

	out_path = "./treated/"

	fig_path= "./figures/"

	rws_save_file = "water_contents_external_calibration.csv" # this will save the output values

	rws_save_fig = "" # not used in the external mode

	paths = (in_liste,in_path,out_path,fig_path,rws_save_file,rws_save_fig)

	switches = ("external","no","no","no")

	rameau(paths,switches,input_properties = ('\t',0))

-----------------------------------
Note on the input file liste
-----------------------------------

The great news about RamEau in Julia is that you can work your file liste in Excel, as it is now a CSV file. It makes it much more pleasant to use, and readable.

If using the "internal" mode, this file liste MUST contain:

	column 1: the file name and extensions, e.g. myspectrum.txt;

	column 2: the name of your product;

	column 3: the water content, if known. If unknow, put 0.0;

	column 4: the spline coefficient for the silicate part. Note: this value is used in the single baseline procedure for the whole spectrum;

	column 5: the spline coefficient for the water part, in case you use the experimental mode with the double baseline fitting procedure (experimental? = "yes" + temperature_laser_correction? = "yes");

	columns 6 to end: the beginning and ends of the BIRs, paired. Please keep the same number of BIRs for all the spectra in one batch.

If using the "external" mode, this file liste MUST contain:

	column 1: the file name and extensions of the references, e.g. myreference.txt;

	column 2: the name of your references;

	column 3: the water content of the references, in wt%;

	column 4: the density of the references, in kg m-3;

	column 5: the file name and extensions of the samples, e.g. mysample.txt;

	column 6: the name of your samples;

	column 7: the estimated density of your samples, in kg m-3.

WARNING: BE SURE THAT THE NUMBER YOU PROVIDE ARE FLOAT NUMBER!

-----------------------------------------------------------------------------
Note on the temperature and excitation line effects corrections
-----------------------------------------------------------------------------

The "internal" mode uses the "long" mode of the tlcorrection function, whereas the "external" mode uses the "hehlen", which takes into account the sample density (see tlcorrection function documentation). This allows to intrisically correct the intensity from density effects.

-----------------------------------------
Note on the experimental mode
-----------------------------------------

The experimental mode contains code for solutions that are currently under development. You may prefer not using it.

However, an interesting feature is provided there, the "double" mode:

When setting the switch experimental? to "double" and combining it with the switch tlcorrection "yes", it allows you to use different smoothing coefficients for the silicate and water signals. In order to use it, you must set the wavenumber of the first ROI for the water band above 2500 cm-1, and the last fo the silicate band below 1600 cm-1 (see the example file for instance). The two different smoothing coefficients are indicated in the dataliste csv file.

-----------------------------------------------------------------------
Note on the use of KRregression baseline fitting instead of GCV splines
-----------------------------------------------------------------------

This is to be used with the internal calibration mode.

Back in 2012 we mostly used the Generalized Cross-Validated splines for fitting the spectral background. However, recent developments show that KRregression or SVMregression may provid better results with less headache for the user (not need to tune the spline coefficient parameter). From experience, using a spline carefully adjusted provides better result. However, using KRregression may provide good results without headache to adjust any parameter. For now this is an experimental feature.

Updates Spetember 2016: A well-adjusted gcvspline usually outperforms the KRregression mode. I advise sticking with the gcvspline for now.
