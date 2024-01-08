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

using Dierckx
using LinearAlgebra
using LsqFit
using Polynomials
using PyCall
using Random
using SparseArrays
using StatsBase
using Statistics

# For PyCall modules
const rampy = PyNULL()

function __init__()
	copy!(rampy, pyimport("rampy"))
end

include("integrale.jl")
include("functions.jl")
include("baseline.jl")
include("bootstrap.jl")
include("tlcorrection.jl")
include("deprecated.jl")
include("ml.jl")
include("peakmeasurement.jl")

#From integrale.jl
export trapz, bandarea

#From functions.jl
export poly, gaussiennes, lorentziennes, pseudovoigts, pearson7, normal_dist
export xshift_inversion, xshift_direct,xshift_correction
export smooth, flipsp, resample, centroid, normalise, despiking

#From baseline.jl
export baseline

#From bootstrap
export bootsample, bootperf

#From tlcorrection
export tlcorrection

#From ml_regressor
export mlregressor, mlexplorer

# From peakmeasurement
export peakmeas

end # module
