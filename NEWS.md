# Spectra.jl News

Copyright (c) 2016 Dr. Charles Le Losq

email: charles.lelosq@anu.edu.au

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.53940.svg)](http://dx.doi.org/10.5281/zenodo.53940)

Licence MIT: see LICENCE.md

As Spectra.jl starts to grow, I will summarize changes in this file starting at version 0.2.0

# Under construction: 0.2.1

Corrections:

- Corrections for removing the warning messages going with Julia 0.5.0

# Stable: 0.2.0

Additions:

- A Ctxremoval function for removing crystal signals from Raman spectra of glasses is added;

- A mlregressor function adding functionalities of SciKit Learn into Spectra.jl for using machine learning in the treatment of spectra has been added.

Improvements:

- Improved Baseline function. IMPORTANT NOTE: the gcvspline and Dspline smoothing coefficients in your software will need to be revised!

- Improved RamEau function with improved internal calibration and critical corrections for the external calibration protocol.

Corrections:

- RamEau external calibration was not converting the water contents from mol/L to wt%. This is now corrected.