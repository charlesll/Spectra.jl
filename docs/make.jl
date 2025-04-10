push!(LOAD_PATH,"../src/")
using Documenter, Spectra

makedocs(sitename="Spectra documentation", 
pages = [
        "index.md",
        "PreProcessing" => "PreProcessing.md",
        "Baseline" => "Baselines.md",
        "Smoothing" => "Smoothing.md",
        "General Functions" => "GeneralFunctions.md",
        "Measurements" => "Measurements.md",
        "Peak Fitting" => "PeakFitting.md",
        "Machine Learning" => "MachineLearning.md",
        "Tips" => "Tips.md",
        "References" => "References.md",
    ],
format = Documenter.HTML(prettyurls = false))

#deploydocs(
#	repo = "github.com/charlesll/Spectra.jl.git",
#	)
#     deps = Deps.pip("mkdocs", "mkdocs-material"),
#     repo = "github.com/charlesll/Spectra.jl.git",
# 	julia  = "1.0",
# 	osname = "linux",
# )
