module MatrixDepot
using GZip, Printf, DelimitedFiles
using LinearAlgebra, SparseArrays, SuiteSparse
import Base: show

export matrixdepot
export info, load, mdopen, mdclose, matrix, rhs, solution, list, metareader, overview
export colval, mtranspose, madjoint
export aliasname, @addgroup, @rmgroup, @modifygroup

include("types.jl")         # common data type definitions
include("higham.jl")          # test matrices
include("regu.jl")               # regularization test problem
include("graph.jl")             # adjacency matrices for graphs
include("data.jl")          # global varaibles and matrix data
include("common.jl")        # main functions
include("download.jl")      # download data from the UF and MM sparse matrix collection
include("datareader.jl")    # read matrix data from local storage
include("matrixmarket.jl")  # read matrix data from local storage

function init()
    GROUP = "group.jl"
    GENERATOR = "generator.jl"

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
            include(joinpath(MY_DEPOT_DIR, file))
        end
    end
    include(joinpath(MY_DEPOT_DIR, GENERATOR))

    downloadindices(MATRIX_DB)
    nothing
end

init()

end # end module
