using Documenter, Spectra

makedocs(sitename="Spectra documentation", format = Documenter.HTML(prettyurls = false))

deploydocs(
	repo = "github.com/charlesll/Spectra.jl.git",
	)
#     deps = Deps.pip("mkdocs", "mkdocs-material"),
#     repo = "github.com/charlesll/Spectra.jl.git",
# 	julia  = "1.0",
# 	osname = "linux",
# )
