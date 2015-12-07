module MatrixDepot
using GZip

import Base: show, search

export matrixdepot, @addgroup, @rmgroup


include("ufreader.jl")
include("common.jl")        # main functions
include("higham.jl")        # test matrices
include("regu.jl")          # regularization test problem
include("data.jl")          # matrix data
include("download.jl")      # download data from the UF sparse matrix collection

const MY_DEPOT_DIR = joinpath(dirname(@__FILE__), "..", "myMatrixDepot")


if !isdir(MY_DEPOT_DIR)
    mkdir(MY_DEPOT_DIR)
    open(string(MY_DEPOT_DIR, "/group.jl"), "w") do f
        write(f, "usermatrixclass = \n Dict( \n \n \n );")
    end
    open(string(MY_DEPOT_DIR, "/generator.jl"), "w") do f
        write(f, "# include your matrix generators below ")
    end
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

end # end module
