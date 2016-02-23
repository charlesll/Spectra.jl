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

function baseline(x::Array{Float64},y::Array{Float64},roi::Array{Float64},basetype::String,p::Array{Float64})
    # First we grab the good roi
    for i = 1:size(roi)[1]
        if i == 1
            interest_index = find(roi[i,1] .<= x[:,1] .<= roi[i,2])
        else
            interest_index = vcat(interest_index,  find(roi[i,1] .<= x[:,1] .<= roi[i,2]))
        end
    end
    interest_x = x[interest_index,1]
    interest_y = y[interest_index,1]
    
    if basetype == "constant"
        return y[:] - minimum(interest_y)
    elseif basetype == "poly"
            return y - poly(p,x) # To be continued... Fitting procedure to add there.
    end
end

    
        