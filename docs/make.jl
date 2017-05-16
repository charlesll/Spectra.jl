using Documenter, Spectra

makedocs()

deploydocs(
    repo = "github.com/charlesll/Spectra.jl.git", 
	julia  = "0.5",
	osname = "linux",
)