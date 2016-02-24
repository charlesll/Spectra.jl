module SpectraJu

include("diffusion.jl")
include("integrale.jl")
include("functions.jl")
include("baseline.jl")

export peak_diffusion, model, IRdataprep, trapz

#From functions.jl
export poly, gaussiennes, normal_dist

export baseline

end # module
