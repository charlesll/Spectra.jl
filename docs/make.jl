using Documenter, Spectra

makedocs(
    modules = [Spectra],
    format = :html,
    sitename = "Spectra.jl",
    pages = Any[
        "Home" => "index.md",
		"Installation" => "Installation.md",
		"Pre-Processing" => "PreProcessing.md",
		"General Functions" => "GeneralFunctions.md",
		"Rameau" => "Rameau.md",
		"Machine Learning Regression" => "MLregressor.md",
		"Peak Fitting" => "PeakFitting.md",
		"Tutorial" => "Tutorial.md",
		"Tips" => "Tips.md",
		"To do" => "ToDo.md",
		"References" => "References.md"
    ],
    doctest = false
)

deploydocs(
    repo = "github.com/charlesll/Spectra.jl.git",
	target = "build",
	julia  = "0.5",
	osname = "linux",
)