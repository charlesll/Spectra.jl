#############################################################################
#Copyright (c) 2016-2019 Charles Le Losq
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

VERSION >= v"1.0.0" && __precompile__()

module Spectra

using DataInterpolations
using Dierckx
using Distributions
using DSP
using ForwardDiff
using LinearAlgebra
using LsqFit
using Measurements
using Optim
using Peaks
using Plots
using Polynomials
using PyCall
using QHull
using Random
using RegularizationTools
using SavitzkyGolay
using SparseArrays
using SpecialFunctions
using Statistics
using StatsBase

# For PyCall modules
const rampy = PyNULL()

function __init__()
    copy!(rampy, pyimport("rampy"))
end

include("integrale.jl")
include("functions.jl")
include("preprocessing.jl")
include("baseline.jl")
include("smoothing.jl")
include("bootstrap.jl")
include("tlcorrection.jl")
include("deprecated.jl")
include("ml.jl")
include("peakmeasurement.jl")
include("fitting.jl")

#From integrale.jl
export trapz

#From functions.jl
export poly, normal_dist
export gaussian, lorentzian, pseudovoigt, pearson7, create_peaks
export funexp, funlog
export invcm_to_nm, nm_to_invcm

# From preprocessing.jl
export xshift_direct, correct_xshift
export flipsp, resample, normalise, despiking, get_portion_interest, extract_signal

#From baseline.jl
export baseline
export arPLS_baseline, drPLS_baseline, als_baseline, rubberband_baseline

# From smoothing.jl
export whittaker, ddmat, smooth

#From bootstrap
export bootsample

#From tlcorrection
export tlcorrection

#From ml_regressor
export mlregressor, mlexplorer

# From peakmeasurement
export peakmeas, centroid, find_peaks, area_peaks

# From fitting
export FitContext, 
    FitResult,
    prepare_context, 
    fit_peaks, 
    print_params, 
    plot_fit, 
    get_peak_results, 
    fit_qNewton, 
    fit_Optim,
    bootstrap

end # module
