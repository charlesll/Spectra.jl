# Tips

In this section are listed various tips for the use of Julia and Spectra:

## Installation

If you see errors messages linked to PyCall, you may have a problem with your environment variable. To solve it, type the following commands in the Julia prompt:

	```julia-repl
	julia> ENV["PYTHON"]=""
	julia> Using Pkg
	julia> Pkg.build("PyCall")
	```

At this point it should work. If yes, you now can type ']' in the Julia repl and then enter:
	```julia-repl
	pkg> add Spectra
	```

## Maintenance

The Julia package ecosystem is frequently evolving. Because of that, it is recommended to update frequently.
Type ']' in the Julia repl and then enter:
	```julia-repl
	pkg> update
	```

## Running Spectra

Always be careful to variable types. Functions will return an error message if you do not enter the good type.
