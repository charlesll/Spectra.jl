#############################################################################
#Copyright (c) 2016 Charles Le Losq
#
#The MIT License (MIT)
#
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the #Software without restriction, including without limitation the rights to use, copy, #modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, #and to permit persons to whom the Software is furnished to do so, subject to the #following conditions:
#
#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, #INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
# baseline.jl contains several "baseline" functions. It is directly dependent on functions.jl.
#
#############################################################################

using JuMP
using Ipopt
using Dierckx

function baseline(x::Array{Float64},y::Array{Float64},roi::Array{Float64},basetype::AbstractString,p::Array{Float64})
    # First we grab the good roi
    interest_index::Array{Int64} = find(roi[1,1] .<= x[:,1] .<= roi[1,2])
    if size(roi)[1] > 1
        for i = 2:size(roi)[1]
            interest_index = vcat(interest_index,  find(roi[i,1] .<= x[:,1] .<= roi[i,2]))
        end
    end
    interest_x = x[interest_index,1]
    interest_y = y[interest_index,1]
    
    if basetype == "poly"
        # The model for fitting baseline to roi signal
        mod = Model(solver=IpoptSolver(print_level=0))
        n::Int = size(p)[1] # number of coefficients
        m::Int = size(interest_x)[1] # number of data

        @variable(mod,p_val[i=1:n])
        setvalue(p_val[i=1:n], p[i])
        @NLobjective(mod,Min,sum{( sum{p_val[i]*interest_x[j]^(i-1), i=1:n} - interest_y[j])^2, j=1:m})
        status = solve(mod)
        println("Solver status: ", status)
        best_p::Vector{Float64} = getvalue(p_val)
        y_calc::Array{Float64} = poly(best_p,x)
        return y[:,1] - y_calc, y_calc
	elseif basetype == "spline"
		spl = Spline1D(interest_x,interest_y,s=p[1])
		y_calc = evaluate(spl,x[:,1])
		return y[:,1] - y_calc, y_calc# To be continued... Fitting procedure to add there.
	else
        error("Not implemented, choose between poly and [to come]")
    end
end

    
        