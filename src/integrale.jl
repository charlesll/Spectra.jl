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
# This files contains functions to perform calculations of numeric integrale values on real dataset x y.
#
#
#############################################################################
"""
	trapz(x::Vector{Tx}, y::Vector{Ty}) where {Tx <: Number, Ty <: Number}

Calculate the area under a curve defined by `x` and `y` values using trapezoidal integration. 
	
# Arguments
- `x::Vector{Tx}`: The x-axis values (e.g., time points or wavelengths), where `Tx <: Number`.
- `y::Vector{Ty}`: The y-axis values (e.g., signal intensities), where `Ty <: Number`.

# Returns
- `area::Tx + Ty`: The computed integral value using the trapezoidal rule.

# Examples

## Example 1: Compute the area under a curve:
```julia
x = [0.0, 1.0, 2.0, 3.0]
y = [0.0, 1.0, 4.0, 9.0]
area = trapz(x, y)
```

## Example 2: Handle edge cases:
```julia
x = [1.0]
y = [2.0]
area = trapz(x, y) # Returns zero since integration requires at least two points.
```

# Notes
- The input vectors `x` and `y` must have the same length; otherwise, an error is raised.
- If the input contains only one point (`length(x) == 1`), the function returns zero as no integration can be performed.

"""
function trapz(x::Vector{Tx}, y::Vector{Ty}) where {Tx<:Number,Ty<:Number}
    # Trapezoidal integration rule
    local n = length(x)
    if (length(y) != n)
        error("Vectors 'x', 'y' must be of same length")
    end
    r = zero(zero(Tx) + zero(Ty))
    if n == 1
        return r
    end
    for i in 2:n
        r += (x[i] - x[i - 1]) * (y[i] + y[i - 1])
    end
    trapz_int = r/2
    return trapz_int
end
