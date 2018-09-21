# The implementations are inspired by MatrixMarket.jl
# https://github.com/JuliaSparse/MatrixMarket.jl
# The MatrixMarket.jl package is licensed under the MIT Expat License:
# Copyright (c) 2013: Viral B. Shah.

# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:

# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
# LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
# OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
# WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

"""
    MatrixDepot

Give access to a wealth of sample and test matrices and accompanying data.
A set of matrices is generated locally (with arguments controlling the special case).
Another set is loaded from one of the publicly accessible matrix collections
`SuiteSparse Matrix Collection` (formerly `University of Florida Matrix Collection`)
and the `Matrix Market Collection`.

Access is like

    using MatrixDepot

    A = matrixdepot("hilb", 10) # locally generated hilbert matrix dimensions (10,10)
    
    A = ("HB/1138_bus")     # named matrix of the SuiteSpares Collection

   or

    md = mdopen("*/bfly")   # named matrix with some extra data
    A = md.A
    co = md.coord
    tx = md("Gname_10.txt")

   or also

    md = mdopen("gravity", 10, false) # localy generated example with rhs and solution
    A = md.A
    b = md.b
    x = md.x
"""
module MatrixDepot
using LinearAlgebra, SparseArrays, SuiteSparse, GZip, Serialization
import Base: show

export matrixdepot
export listnames, listdir, listdata, listgroups, mdlist, mdinfo, metasymbols, mdopen
export @addgroup, @rmgroup, @modifygroup
#export colval, mtranspose, madjoint
# further exports (for predicate functions) in `logical.jl`

include("types.jl")         # common data type definitions
include("higham.jl")        # test matrices
include("regu.jl")          # regularization test problem
include("graph.jl")         # adjacency matrices for graphs
include("data.jl")          # global variables and matrix data
include("common.jl")        # main functions
include("logical.jl")       # operations on patterns and predicates
include("download.jl")      # download data from the UF and MM sparse matrix collection
include("datareader.jl")    # read matrix data from local storage
include("matrixmarket.jl")  # read matrix data from local storage
include("markdown.jl")      # construct MD objects

function init(;ignoredb::Bool=false)
    GROUP = "group.jl"
    GENERATOR = "generator.jl"

    if !isdir(DATA_DIR)
        mkpath(DATA_DIR)
    end

    if !isdir(MY_DEPOT_DIR)
        mkpath(MY_DEPOT_DIR)
        open(joinpath(MY_DEPOT_DIR, GROUP), "w") do f
            write(f, "usermatrixclass = Dict(\n);")
        end
        open(joinpath(MY_DEPOT_DIR, GENERATOR), "w") do f
            write(f, "# include your matrix generators below \n")
        end
        println("created dir $MY_DEPOT_DIR")
    end
    
    for file in readdir(MY_DEPOT_DIR)
        if endswith(file, ".jl") && file != GENERATOR
            println("include $file for user defined matrix generators")
            include(joinpath(MY_DEPOT_DIR, file))
        end
    end
    include(joinpath(MY_DEPOT_DIR, GENERATOR))
    println("verify download of index files...")
    downloadindices(MATRIX_DB, ignoredb=ignoredb)
    println("used remote site is $(uf_remote.params.indexurl)")
    println("populating internal database...")
    nothing
end

# will be called automatically once after `using`, `import`, `require`.
function __init__()
    try init() catch ex; @warn "exception during initialization: '$ex'"
    end
end

end # end module
