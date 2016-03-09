# Part of the code of the function denseread come from 
#The MatrixMarket.jl package which is licensed under the MIT Expat License:

#Copyright (c) 2013: Viral B. Shah.

#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

function denseread(filename::AbstractString)
    file = open(filename, "r")
    firstline = chomp(readline(file))
    tokens = split(firstline)
    if tokens[1] != "%%MatrixMarket"
        throw(ParseError(string("Not a valid MatrixMarket header:", ll)))
    end
    (head1, rep, field, symm) = map(lowercase, tokens[2:5])
    eltype = field == "real" ? Float64 :
                  field == "complex" ? Complex128 :
                  field == "pattern" ? Bool :
                  throw(ParseError("Unsupported field $field (only real and complex are supported)"))
    # Skip all comments and empty lines
    ll   = readline(file)
    while length(chomp(ll))==0 || (length(ll) > 0 && ll[1] == '%')
        ll = readline(file)
    end
    return readdlm(file, eltype)
end
