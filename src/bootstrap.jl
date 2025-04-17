#############################################################################
#Copyright (c) 2016-2025 Charles Le Losq
#
#The MIT License (MIT)
#
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the #Software without restriction, including without limitation the rights to use, copy, #modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, #and to permit persons to whom the Software is furnished to do so, subject to the #following conditions:
#
#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, #INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
# functions.jl contains several mathematic functions
#
#
#############################################################################

"""
	bootsample(x, y; boottype::String="np", ese=nothing)

# Inputs
- `x`: the x axis. It can have multiple columns.
- `y`: the y axis. It can have multiple columns.

# Options
- `boottype::String = "np"`: type of bootstrapping
    - "np": non-parametric bootstrapping. Data resampled with replacement. 
    - "p": parametric bootstrapping. Data resampled from the Normal(y, ese) distribution.
- `ese = nothing`: standard errors on `y`.

# Outputs
- `b_x_f`: bootstrapped x sample
- `b_y_f`: bootstrapped y sample

The bootstrap function can be embedded in a for loop, and will each time produce a different dataset. Performing K times the bootstrapping and fitting each time the model will allow to estimate the error distribution on the peak parameters. This technic has the advantage of making no prior assumption on the probability distribution functions of parameters errors. However, it is more time consuming that using the covariance matrix.

# References
- Efron, B. 1979. “Bootstrap Methods: Another Look at the Jackknife.” The Annals of Statistics 7 (1): 1–26.
- Efron, Bradley. 1981. “Nonparametric Estimates of Standard Error: The Jackknife, the Bootstrap and Other Methods.” Biometrika 68 (3): 589–99. doi:10.1093/biomet/68.3.589.
- Efron, B., and Tibshirani, R. 1994. An Introduction to the Bootstrap. CRC press.

"""
function bootsample(x, y; boottype::String="np", ese=nothing)
    testx = size(x, 1)
    testy = size(y, 1)
    if ese == nothing && boottype == "p"
        error("You must provide standard errors for parametric bootstrapping.")
    elseif ese != nothing
        teste = size(ese, 1)
        if teste != testx
            error("Check the size of x, y and ese.")
        end
    end

    if boottype == "np" && testx == testy
        vect = collect(1:size(x)[1]) # for real bootstrapping
        idx = sample(vect, size(vect)[1]; replace=true) #resampling data with replacement...
        b_x_f = x[idx, :] # we pick the right x
        b_y_f = y[idx, :] # we pick the right y
        if ese != nothing 
            # if the user provides errors, we also return a sample of them
            b_e_f = ese[idx, :]
            return b_x_f, b_y_f, b_e_f
        else
            return b_x_f, b_y_f
        end
    elseif boottype == "p" && testx == testy && testx == teste
        b_x_f = x
        b_y_f = ones(size(y))
        b_y_f = randn!(b_y_f[:, :]) .* ese[:, :] + y[:, :]
        return b_x_f, b_y_f
    else
        error(
            "Something is wrong. Check size of x, y and ese as well as the boottype (either p or np). Providing an ese array is mandatory for parametric bootstrapping",
        )
    end
end