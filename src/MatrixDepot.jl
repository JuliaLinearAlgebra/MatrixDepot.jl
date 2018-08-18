module MatrixDepot
using GZip, Printf, DelimitedFiles
using LinearAlgebra, SparseArrays, SuiteSparse

import Base: show

export matrixdepot, @addgroup, @rmgroupa, info, load, mdread


include("types.jl")         # common data type definitions
include("higham.jl")          # test matrices
include("regu.jl")               # regularization test problem
include("graph.jl")             # adjacency matrices for graphs
include("data.jl")          # global varaibles and matrix data
include("common.jl")        # main functions
include("download.jl")      # download data from the UF and MM sparse matrix collection
include("datareader.jl")    # read matrix data from local storage
include("matrixmarket.jl")  # special handling for MM


const MY_DEPOT_DIR = joinpath(dirname(@__FILE__), "..", "myMatrixDepot")

function init()

    if !isdir(MY_DEPOT_DIR)
        mkdir(MY_DEPOT_DIR)
        open(string(MY_DEPOT_DIR, "/group.jl"), "w") do f
            write(f, "usermatrixclass = \n Dict( \n \n \n );")
        end
        open(string(MY_DEPOT_DIR, "/generator.jl"), "w") do f
            write(f, "# include your matrix generators below ")
        end
        println("created dir $MY_DEPOT_DIR")
    end

    files = Set(readdir(MY_DEPOT_DIR))
    delete!(files, "generator.jl")
    if isdir(MY_DEPOT_DIR)
        for file in files
            if split(file, '.')[2] == "jl"
                include("$(MY_DEPOT_DIR)/$(file)")
            end
        end
        include(string(MY_DEPOT_DIR, "/generator.jl"))
    end

    downloadindices(MATRIX_DB)
    nothing
end

init()

end # end module
