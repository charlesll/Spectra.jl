push!(LOAD_PATH, "../src/")
using Documenter, DocumenterTools, Spectra, Literate, Plots#, Turing

# Code for examples from https://github.com/oxfordcontrol/COSMO.jl/blob/master/docs/make.jl

@info "Building example problems..."

# utility function from https://github.com/JuliaOpt/Convex.jl/blob/master/docs/make.jl
fix_math_md(content) = replace(content, r"\$\$(.*?)\$\$"s => s"```math\1```")
fix_suffix(filename) = replace(filename, ".jl" => ".md")
function postprocess(cont)
    """
    The source files for all examples can be found in [/examples](https://github.com/charlesll/Spectra.jl/tree/master/examples/).
    """ * cont
end

# find all example source files
exclude_files = [];
example_path = joinpath(@__DIR__, "../examples/")
build_path = joinpath(@__DIR__, "src", "examples/")
files = readdir(example_path)
filter!(x -> endswith(x, ".jl"), files)
filter!(x -> !in(x, exclude_files), files)

for file in files
    Literate.markdown(
        example_path * file,
        build_path;
        preprocess=fix_math_md,
        postprocess=postprocess,
        documenter=true,
        credit=true,
    )
end

examples_nav = fix_suffix.("./examples/" .* files)

@info "Makeing documentation..."

makedocs(;
    sitename="Spectra documentation",
    authors="Charles Le Losq and contributors.",
    pages=[
        "index.md",
        "Data Processing" => "PreProcessing.md",
        "Measurements" => "Measurements.md",
        "Peak Fitting" => "PeakFitting.md",
        "Machine Learning" => "MachineLearning.md",
        "Helper Functions" => "HelperFunctions.md",
        "Tutorials" => examples_nav,
        "Tips" => "Tips.md",
        "References" => "References.md",
    ],
    format=Documenter.HTML(; prettyurls=false),
)

#deploydocs(
#	repo = "github.com/charlesll/Spectra.jl.git",
#	)
#     deps = Deps.pip("mkdocs", "mkdocs-material"),
#     repo = "github.com/charlesll/Spectra.jl.git",
# 	julia  = "1.0",
# 	osname = "linux",
# )
