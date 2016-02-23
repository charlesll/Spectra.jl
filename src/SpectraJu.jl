module SpectraJu

include("diffusion.jl")
include("integrale.jl")
include("functions.jl")
include("baseline.jl")

export peak_diffusion, model, IRdataprep, trapz

export poly

export baseline

end # module
