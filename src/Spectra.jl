module Spectra

include("diffusion.jl")
include("integrale.jl")
include("functions.jl")
include("baseline.jl")

#From integrale.jl
export trapz

#From diffusion.jl
export peak_diffusion, model, IRdataprep 

#From functions.jl
export poly, gaussiennes, normal_dist, bootstrap

#From baseline.jl
export baseline

end # module
