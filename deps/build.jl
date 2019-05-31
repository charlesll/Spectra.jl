using PyCall, Conda

println("Running build.jl for the Spectra package.")

# Change that to whatever packages you need.
const PACKAGES = ["rampy>=0.4.4"]

pyimport_conda("pip", "pip")

#Conda.add("pip")
run(`$(PyCall.python) -m pip install $(PACKAGES)`)#
#run(`pip install --upgrade gcvspline`)

# Change that to whatever packages you need.
# const PACKAGES = ["gcvspline"]

#= # Use eventual proxy info
proxy_arg=String[]
if haskey(ENV, "http_proxy")
    push!(proxy_arg, "--proxy")
    push!(proxy_arg, ENV["http_proxy"])
end

# Import pip
try
    @pyimport pip
catch
    # If it is not found, install it
    println("Pip not found on your system. Downloading it.")
    get_pip = joinpath(dirname(@__FILE__), "get-pip.py")
    download("https://bootstrap.pypa.io/get-pip.py", get_pip)
    run(`$(PyCall.python) $(proxy_arg) $get_pip --user`)
end

println("Installing required python packages using pip")
run(`$(PyCall.python) $(proxy_arg) -m conda install --upgrade pip setuptools`)
run(`$(PyCall.python) $(proxy_arg) -m pip install $(PACKAGES)`) =#
