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

VERSION >= v"0.4.0" && __precompile__()

module Spectra

using StatsBase
using PyPlot

# some initial setup for calling the GCVSPL.f library
unixpath = "../deps/src/gcvspline/libgcvspl"
winpath = "../deps/bin$WORD_SIZE/libgcvspl" # let it there as an example but I did not tried yet any build on Windows... TODO
const gcvspl = joinpath(dirname(@__FILE__), @unix? unixpath : winpath)

function __init__()
    # Ensure library is available.
    if (Libdl.dlopen_e(gcvspl) == C_NULL)
        error("GCVSPL not properly installed. Run Pkg.build(\"Spectra\"). Windows auto-build is not setup, you might want to build the library manually.")
    end
end

include("diffusion.jl")
include("integrale.jl")
include("functions.jl")
include("baseline.jl")
include("bootstrap.jl")
include("long.jl")

#From integrale.jl
export trapz, gaussianarea

#From diffusion.jl
export peak_diffusion, model, IRdataprep

#From functions.jl
export poly, gaussiennes, lorentziennes, pseudovoigts, pearson7, normal_dist

#From baseline.jl
export baseline

#From gcvspl.jl
export gcvspl, splderivative

#From bootstrap
export bootsample, bootperf

#From long
export long

end # module
