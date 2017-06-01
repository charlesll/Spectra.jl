# Spectra.jl News

Copyright (c) 2016-2017 Dr. Charles Le Losq

email: charles.lelosq@anu.edu.au

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.53940.svg)](http://dx.doi.org/10.5281/zenodo.53940)

Licence MIT: see LICENCE.md

As Spectra.jl starts to grow, I will summarize changes in this file starting at version 0.2.0

# 0.3.4

- Wrapping of the GCVSPL.f library is dropped, we now call gcvspline, the Python package that wraps GCVSPL.f. This now should solve any trouble with Windows users! No change of API, except the drop of the gcvspl_julia and splder_julia functions. This should be transparent for users.

- smooth function added, see doc; smoothing with gcvspline is provided.

- Numerous tests added.

- Various examples added.

- Spectra works with Julia v0.6: test on Travis OK.

# 0.3.3


- Critical change in peakhw that was returning the full width of the peaks instead of the half-width. This is now corrected. Name of peakhw is changed in peakmeas. This also outputs the peak intensity and centroid.

- Merged pull request #6 from staticfloat/updated_ci_url

# 0.3.2

- Various slight corrections in the code

- Documentation is now generated using Documenter.jl

# 0.3.1

Corrections:

- In the documentation of the rameau function, the prediction coefficient was badly indicated. This is an array containing the coefficient and its error, not a single float value.

- /deps/build.jl was badly written. This is now corrected.

# 0.3.0

Critical changes in requirements:

- Compat removed

- Julia >= 0.5 (dropping 0.4.x)

Additions:

- mmap_switch option in rameau added to overpass the problem of memory mapping error when calling readcsv/readdlm in virtual environments

- roi_hf_external is a new option that allows to setup where the baseline is fitted when using the "External" calibration methods

- scaler option in mlregressor allows to choose between two ways of standardizing the data for machine learning

Corrections:

- A few corrections for better support of Julia versions >= 0.5

- A few corrections and documentation improvements

# 0.2.3 (includes 0.2.2; stable release)

Additions:

- Addition of the double baseline feature in Rameau in the experimental mode

Corrections:

- Corrections for removing the warning messages going with Julia 0.5.0

- Use of Compat and Compat.String to avoid warning in Julia 0.4


# 0.2.1

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
