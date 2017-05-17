using Documenter, Spectra

makedocs()

deploydocs(
    deps = Deps.pip("mkdocs", "mkdocs-material"),
    repo = "github.com/charlesll/Spectra.jl.git", 
	julia  = "0.5",
	osname = "linux",
)