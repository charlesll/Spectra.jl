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

VERSION >= v"0.5.0" && __precompile__()

module Spectra

using StatsBase
using PyPlot
using LsqFit
using PyCall
using Dierckx
using Ipopt
using JuMP
using NMF

# For PyCall modules
const preprocessing = PyNULL()
const model_selection = PyNULL()
const decomposition = PyNULL()
const kernel_ridge = PyNULL()
const svm = PyNULL()
const linear_model = PyNULL()
const gaussian_process = PyNULL()
const pygcvspl = PyNULL()

function __init__()
	
	copy!(preprocessing, pyimport_conda("sklearn.preprocessing", "scikit-learn"))
	copy!(model_selection, pyimport_conda("sklearn.model_selection", "scikit-learn"))
	copy!(decomposition, pyimport_conda("sklearn.decomposition", "scikit-learn"))
	copy!(kernel_ridge, pyimport_conda("sklearn.kernel_ridge", "scikit-learn"))
	copy!(svm, pyimport_conda("sklearn.svm","scikit-learn"))
	copy!(gaussian_process, pyimport_conda("sklearn.gaussian_process", "scikit-learn"))
	copy!(linear_model, pyimport_conda("sklearn.linear_model", "scikit-learn"))
	copy!(pygcvspl, pyimport("gcvspline"))

end

include("diffusion.jl")
include("integrale.jl")
include("functions.jl")
include("baseline.jl")
include("bootstrap.jl")
include("tlcorrection.jl")
include("rameau.jl")
include("deprecated.jl")
include("ml_regressor.jl")
include("ctxremoval.jl")
include("peakmeasurement.jl")

#From integrale.jl
export trapz, bandarea

#From diffusion.jl
export peak_diffusion, model, IRdataprep

#From functions.jl
export poly, polyfit, gaussiennes, lorentziennes, pseudovoigts, pearson7, normal_dist, xshift_inversion, xshift_direct,xshift_correction, SavitzkyGolayFilter, smooth

#From baseline.jl
export baseline, whitsmdd

#From bootstrap
export bootsample, bootperf

#From tlcorrection
export tlcorrection

#From rameau
export rameau

#From ml_regressor
export mlregressor

#From ctxremoval
export ctxremoval

# From peakmeasurement
export peakmeas

end # module
