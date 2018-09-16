module MatrixDepot
using GZip, Printf, DelimitedFiles
using LinearAlgebra, SparseArrays, SuiteSparse
import Base: show

export matrixdepot
export metareader
export mdlist, listdir, listdata, mdinfo, metadata, mdopen
export @addgroup, @rmgroup, @modifygroup
#export colval, mtranspose, madjoint
# further exports (for predicate functions) in `logical.jl`

include("types.jl")         # common data type definitions
include("higham.jl")        # test matrices
include("regu.jl")          # regularization test problem
include("graph.jl")         # adjacency matrices for graphs
include("data.jl")          # global varaibles and matrix data
include("common.jl")        # main functions
include("logical.jl")       # operations on patterns and predicates
include("download.jl")      # download data from the UF and MM sparse matrix collection
include("datareader.jl")    # read matrix data from local storage
include("matrixmarket.jl")  # read matrix data from local storage
include("markdown.jl")      # construct MD objects

function init()
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
    downloadindices(MATRIX_DB)
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
