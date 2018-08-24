# Part of the code of the function denseread come from 
#The MatrixMarket.jl package which is licensed under the MIT Expat License:

#Copyright (c) 2013: Viral B. Shah.

#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

function mmread(filename::AbstractString)
    open(filename, "r") do file
        line = lowercase(readline(file))
        #println("tokens: ", line)
        tokens = split(line)
        if tokens[1] != "%%matrixmarket"
            throw(Meta.ParseError(string("Not a valid MatrixMarket header:", line)))
        end
        line = readline(file)
        while length(line) == 0 || line[1] == '%'
            line = readline(file)
        end
        if tokens[2] == "matrix"
            mmread_matrix(file, line, tokens[3:end]...)
        else
            throw(Meta.ParseError(string("Unsupported type: ", line)))
        end
    end
end

lt(r::Matrix) = (a, b) -> r[2,a] < r[2,b] || ( r[2,a] == r[2,b] && r[1,a] < r[1,b] )

function mmread_matrix(file::IO, line, rep, field, symm)
    FMAP = Dict("real" => Float64, "complex" => ComplexF64,
                "integer" => Int64, "pattern" => Bool)

    T = get(FMAP, field) do
         throw(Meta.ParseError("Unsupported field $field (only real/complex/pattern)"))
    end

    sparse = rep == "coordinate"

    if sparse
        m, n, nz = parseint(line)
        rv = Vector{Int}(undef, nz)
        cv = Vector{Int}(undef, nz)
        vv = Vector{T}(undef, nz)
        for i = 1:nz
            line = readline(file)
            parseline!(rv, cv, vv, i, line)
        end
        result = mksparse(m, n, rv, cv, vv)
    else
        n, m = parseint(line)
        result = zeros(T, m, n)
        for c = 1:n
            for r = 1:m
                line = readline(file)
                x = parseline(line, T)
                result[r,c] = x[1]
            end
        end
    end
                
    if symm == "symmetric"
        Symmetric(result, :L)
    elseif symm == "skew-symmetric"
        result - transpose(result)
    elseif symm == "general"
        result
    elseif symm == "hermitian"
        Hermitian(result, :L)
    else
        throw(Meta.ParseError("Unsupported qualifier $symm (only general/symmetric/skew-symmetric/hermitian)"))
    end
end

function parseline!(rv::Vector{Int}, cv::Vector{Int}, vv::Vector{T},
                    i::Int, line::AbstractString) where T

    tokens = split(line)
    2 <= length(tokens) <= 4 || error("parsing $line for (Int, Int, $T)")
    #println("parseline($line, $args) tokens=$tokens")
    rv[i] = parse(Int, tokens[1])
    cv[i] = parse(Int, tokens[2])
    vv[i] = parsevalue(tokens, T)
end

function parseline(line::AbstractString, T::Type)
    tokens = split(line)
    parsevalue(tokens, T)
end

parsevalue(tokens::Vector{<:AbstractString},T::Type{<:Bool}) = true
parsevalue(tokens::Vector{<:AbstractString},T::Type{<:Real}) = parse(T, tokens[end])
function parsevalue(tokens::Vector{<:AbstractString},T::Type{<:Complex})
    R = real(T)
    parse(R, tokens[end-1]) + parse(R, tokens[end])*im
end

function parseint(line::AbstractString)
    tokens = split(line)
    parse.(Int, tokens)
end

function mksparse(m::Integer, n::Integer, rv::AbstractVector{Ti}, cv::AbstractVector{Ti},
                  vv::AbstractVector{Tv}) where {Ti<:Integer,Tv}

    nz = length(rv)
    p = sortperm(rv)
    q = sortperm(cv[p])
    p = q[p]
    colptr = Vector{Int}(undef, n+1)
    colp = 0
    for i = 1:nz
        while colp < cv[p[i]]
            colp += 1
            colptr[colp] = i
        end
    end
    while colp <= n
        colp += 1
        colptr[colp] = nz + 1
    end
    SparseMatrixCSC{Tv,Ti}(m, n, colptr, rv[p], vv[p])
end
